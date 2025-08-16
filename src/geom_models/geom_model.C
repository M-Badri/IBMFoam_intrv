/*---------------------------------------------------------------------------*\
License 

   IBMFoam is distributed under the GNU Lesser General Public License (LGPL).
   
   You are free to copy and share this license text in its original form. 
   Modifying the wording of the license itself is not permitted.
   
   This license incorporates the rights and obligations of the 
   GNU General Public License (GPL) v3, 
   along with the additional permissions granted under the LGPL terms.
   
   A copy of the GNU Lesser General Public License should have been provided 
   with IBMFoam. If you did not receive one, it can be found online at:
      <http://www.gnu.org/licenses/lgpl.html>

InNamspace
    Foam
\*---------------------------------------------------------------------------*/
#include "geom_model.H"

using namespace Foam;

//---------------------------------------------------------------------------//
geom_model::geom_model
(
    const  fvMesh&   mesh,
    const contactType cType,
    scalar  thrSurf
)
:
contactType_(cType),
mesh_(mesh),
ibPartialVolume_(Pstream::nProcs(), 0),
owner_(0),
cellToStartInCreateIB_(0),
thrSurf_(thrSurf),
intSpan_(2.0),
sdBasedLambda_(false),
curMeshBounds_(mesh_.points(),false),
M_(0.0),
M0_(0.0),
nCells_(0),
CoM_(vector::zero),
I_(symmTensor::zero),
bBox_(std::make_shared<boundBox>()),
dC_(0.0),
rhoS_("rho",dimensionSet(1,-3,0,0,0,0,0),1.0)
{
    surfCells_.setSize(Pstream::nProcs());
    intCells_.setSize(Pstream::nProcs());
}
geom_model::~geom_model()
{
}
//---------------------------------------------------------------------------//
void geom_model::calculateGeometricalProperties
(
    volScalarField& body
)
{
    M_      = scalar(0);
    setCoM();
    I_      = symmTensor::zero;

    nCells_ = label(0);

    addToMAndI(body,surfCells_[Pstream::myProcNo()]);
    addToMAndI(body,intCells_[Pstream::myProcNo()]);

    nCells_ = intCells_[Pstream::myProcNo()].size() + surfCells_[Pstream::myProcNo()].size();


    reduce(M_, sumOp<scalar>());

    reduce(I_,  sumOp<symmTensor>());

    reduce(nCells_, sumOp<label>());

}
//---------------------------------------------------------------------------//
void geom_model::calculateGeometricalPropertiesParallel
(
    volScalarField& body
)
{

    M_      = scalar(0);

    setCoM();
    I_      = symmTensor::zero;

    nCells_ = label(0);

    addToMAndI(body,surfCells_[Pstream::myProcNo()]);
    addToMAndI(body,intCells_[Pstream::myProcNo()]);
    nCells_ = intCells_[Pstream::myProcNo()].size() + surfCells_[Pstream::myProcNo()].size();

}
//---------------------------------------------------------------------------//
void geom_model::addToMAndI
(
    volScalarField& body,
    DynamicLabelList& labelCellLst
)
{
    forAll (labelCellLst,cell)
    {
        label cellI  = labelCellLst[cell];

        scalar Mi    = body[cellI]*rhoS_.value()*mesh_.V()[cellI];

        M_      += Mi;
        scalar xLoc = mesh_.C()[cellI].x() - CoM_.x();
        scalar yLoc = mesh_.C()[cellI].y() - CoM_.y();
        scalar zLoc = mesh_.C()[cellI].z() - CoM_.z();

        // add to I_
        I_.xx() += Mi*(yLoc*yLoc + zLoc*zLoc);
        I_.yy() += Mi*(xLoc*xLoc + zLoc*zLoc);
        I_.zz() += Mi*(xLoc*xLoc + yLoc*yLoc);

        I_.xy() -= Mi*(xLoc*yLoc);
        I_.xz() -= Mi*(xLoc*zLoc);
        I_.yz() -= Mi*(yLoc*zLoc);
    }
}
//---------------------------------------------------------------------------//
void geom_model::compute_body_charPars()
{
    forAll (surfCells_[Pstream::myProcNo()],sCellI)
    {
        label cellI = surfCells_[Pstream::myProcNo()][sCellI];
        dC_ = max(dC_,mag(CoM_-mesh_.C()[cellI]));
    }
    M0_ = M_;
    reduce(dC_, maxOp<scalar>());
}
//---------------------------------------------------------------------------//
void geom_model::resetBody(volScalarField& body)
{
    forAll (intCells_[Pstream::myProcNo()],cellI)
    {
        body[intCells_[Pstream::myProcNo()][cellI]] = 0;
    }
    forAll (surfCells_[Pstream::myProcNo()],cellI)
    {
        body[surfCells_[Pstream::myProcNo()][cellI]] = 0;
    }

    surfCells_[Pstream::myProcNo()].clear();
    intCells_[Pstream::myProcNo()].clear();
    }
//---------------------------------------------------------------------------//
bool geom_model::isBBoxInMesh()
{
    boundBox ibBound(getBounds());

    forAll(geometricD,dir)
    {
        if(geometricD[dir] == 1)
        {
            if(!(curMeshBounds_.max()[dir] >= ibBound.min()[dir]
                && curMeshBounds_.min()[dir] <= ibBound.max()[dir]))
            {
                return false;
            }
        }
    }
    return true;
}
//---------------------------------------------------------------------------//
DynamicList<label> geom_model::getPotentSurfCells
(
    volScalarField& body,
    HashTable<bool, label, Hash<label>>& cellInside,
    List<labelList>& cellPoints
)
{
    const labelList foundCells = cellInside.toc();
    DynamicLabelList potentSurfCells(foundCells.size()/2);
    forAll(foundCells, cellI)
    {
        label cCell = foundCells[cellI];

        if(cellInside[cCell])
        {
            const labelList& neigh = cachedNeighbours_()[cCell];

            bool anyOutside = false;
            forAll(neigh, neighI)
            {
                if(!cellInside[neigh[neighI]])
                {
                    anyOutside = true;
                    break;
                }
            }

            if(anyOutside)
            {
                potentSurfCells.append(cCell);
            }
            else
            {
                intCells_[Pstream::myProcNo()].append(cCell);
                ibPartialVolume_[Pstream::myProcNo()] += 1;
                body[cCell] = 1.0;
            }
        }
        else
        {
            potentSurfCells.append(cCell);
        }
    }
    return potentSurfCells;
}
//---------------------------------------------------------------------------//
void geom_model::correctSurfCells
(
    volScalarField& body,
    DynamicLabelList& potentSurfCells,
    HashTable<bool, label, Hash<label>>& cellInside,
    List<labelList>& cellPoints
)
{
    HashTable<bool, label, Hash<label>> verticesStatus(potentSurfCells.size()*6);
    forAll(potentSurfCells, cellLabel)
    {
        label cCell = potentSurfCells[cellLabel];
        scalar cBody(0);

        if(cellInside[cCell])
        {
            cBody += 0.5;
        }

        const labelList& cVerts = cellPoints[cCell];
        scalar rVInSize(0.5/cVerts.size());
        forAll(cVerts, vertI)
        {
            if(!verticesStatus.found(cVerts[vertI]))
            {
                bool vertexInside = pointInside(mesh_.points()[cVerts[vertI]]);
                verticesStatus.insert
                (
                    cVerts[vertI],
                    vertexInside
                );

                if(vertexInside)
                {
                    cBody += rVInSize;
                }
            }
            else if(verticesStatus[cVerts[vertI]])
            {
                cBody += rVInSize;
            }
        }

        if (cBody > thrSurf_)
        {
            if (cBody > (1.0-thrSurf_))
            {
                intCells_[Pstream::myProcNo()].append(cCell);
            }
            else
            {
                surfCells_[Pstream::myProcNo()].append(cCell);
            }
            ibPartialVolume_[Pstream::myProcNo()] += 1;

            body[cCell] += cBody;


            body[cCell] = min(max(0.0,body[cCell]),1.0);
        }
    }
}
//---------------------------------------------------------------------------//
List<std::shared_ptr<boundBox>> geom_model::getBBoxes()
{
    boundBox cBBox = getBounds();
    bBox_->min() = cBBox.min();
    bBox_->max() = cBBox.max();

    List<std::shared_ptr<boundBox>> retList(1);
    retList[0] = bBox_;
    return retList;
}
//---------------------------------------------------------------------------//
