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
#include "sphere_body.H"

using namespace Foam;

//---------------------------------------------------------------------------//
// create immersed body for convex body
void sphere_body::create_immersed_body
(
    volScalarField& body,
    Field<label>& octreeField,
    List<labelList>& cellPoints
)
{

    intCells_[Pstream::myProcNo()].clear();
    surfCells_[Pstream::myProcNo()].clear();
    ibPartialVolume_[Pstream::myProcNo()] = 0;

    if(!isBBoxInMesh())
    {
        return;
    }

    label cellInIB = getCellInBody(octreeField);
    if(cellInIB == -1)
    {
        return;
    }

    const pointField& cp = mesh_.C();

    autoPtr<DynamicLabelList> nextToCheck(
        new DynamicLabelList(1,cellInIB));
    autoPtr<DynamicLabelList> auxToCheck(
        new DynamicLabelList);

    label tableSize = 128;
    if(cachedNeighbours_.valid())
    {
        tableSize = cachedNeighbours_().toc().size()*1.5;
    }
    else
    {
        cachedNeighbours_ = new HashTable<const labelList&, label, Hash<label>>;
    }

    HashTable<bool, label, Hash<label>> cellInside(tableSize);

    label iterCount(0);label iterMax(mesh_.nCells());
    while (nextToCheck().size() > 0 and iterCount++ < iterMax)
    {
        auxToCheck().clear();
        forAll (nextToCheck(),cellToCheck)
        {
            label cCell = nextToCheck()[cellToCheck];
            if (!cellInside.found(cCell))
            {
                iterCount++;

                if(pointInside(cp[cCell]))
                {
                    cellInside.set(cCell, true);

                    if(cachedNeighbours_.valid() && cachedNeighbours_().found(cCell))
                    {
                        auxToCheck().append(cachedNeighbours_()[cCell]);
                    }
                    else
                    {
                        const labelList& neigh = mesh_.cellCells(cCell);
                        cachedNeighbours_().insert(cCell, neigh);
                        auxToCheck().append(neigh);
                    }
                }
                else
                {
                    cellInside.set(cCell, false);
                }
            }
        }
        const autoPtr<DynamicLabelList> helpPtr(nextToCheck.ptr());
        nextToCheck.set(auxToCheck.ptr());
        auxToCheck = helpPtr;
    }

    DynamicLabelList keyToErase;
    for(auto it = cachedNeighbours_().begin(); it != cachedNeighbours_().end(); ++it)
    {
        if(!cellInside.found(it.key()))
        {
            keyToErase.append(it.key());
        }
    }
    cachedNeighbours_().erase(keyToErase);

    DynamicLabelList potentSurfCells =
        getPotentSurfCells(
            body,
            cellInside,
            cellPoints
        );

    correctSurfCells
    (
        body,
        potentSurfCells,
        cellInside,
        cellPoints
    );

    if(intCells_[Pstream::myProcNo()].size() > 0)
    {
        cellToStartInCreateIB_ = min(intCells_[Pstream::myProcNo()]);
    }
}
//---------------------------------------------------------------------------//
label sphere_body::getCellInBody
(
    Field<label>& octreeField
)
{

    labelHashSet checkedCells;

    const pointField& cp = mesh_.C();

    if(cellToStartInCreateIB_ >= octreeField.size())
        cellToStartInCreateIB_ = 0;

    autoPtr<DynamicLabelList> nextToCheck(
        new DynamicLabelList(1,cellToStartInCreateIB_));
    autoPtr<DynamicLabelList> auxToCheck(
        new DynamicLabelList);

    label iterCount(0);label iterMax(mesh_.nCells());

    while (nextToCheck().size() > 0 and iterCount < iterMax)
    {
        auxToCheck().clear();
        forAll (nextToCheck(),cellToCheck)
        {
            if (!checkedCells.found(nextToCheck()[cellToCheck]))
            {
                checkedCells.insert(nextToCheck()[cellToCheck]);
                iterCount++;

                if(pointInside(cp[nextToCheck()[cellToCheck]]))
                {
                    return nextToCheck()[cellToCheck];
                }
                else
                {
                    auxToCheck().append(mesh_.cellCells()[nextToCheck()[cellToCheck]]);
                }
            }
        }
        const autoPtr<DynamicLabelList> helpPtr(nextToCheck.ptr());
        nextToCheck.set(auxToCheck.ptr());
        auxToCheck = helpPtr;
    }
    return -1;
}
//---------------------------------------------------------------------------//
void sphere_body::synchronPos(label owner)
{
    PstreamBuffers pBufs(Pstream::commsTypes::nonBlocking);

    owner = (owner == -1) ? owner_ : owner;

    if (owner == Pstream::myProcNo())
    {
        for (label proci = 0; proci < Pstream::nProcs(); proci++)
        {
            UOPstream send(proci, pBufs);
            send << position_;
        }
    }

    pBufs.finishedSends();

    UIPstream recv(owner, pBufs);
    vector pos (recv);


    position_ = pos;
}
//---------------------------------------------------------------------------//
boolList sphere_body::pointInside(pointField pointI)
{
    boolList inside(pointI.size());

    forAll(pointI,point)
    {
        inside[point] = mag(position_-pointI[point]) < radius_;
    }

    return inside;
}
//---------------------------------------------------------------------------//
bool sphere_body::pointInside(point pointI)
{
    return mag(position_-pointI) < radius_;
}
//---------------------------------------------------------------------------//
pointField sphere_body::sampleSurfacePoints()
{
    pointField returnField(6);

    vector a(1,0,0);
    vector b(0,1,0);
    vector c(0,0,1);

    List<vector> listV(3);
    listV[0] = a;
    listV[1] = b;
    listV[2] = c;

    forAll (listV,v)
    {
        returnField[v] = position_ + listV[v] * radius_;
    }
    forAll (listV,v)
    {
        returnField[3+v] = position_ - listV[v] * radius_;
    }

    return returnField;
}
