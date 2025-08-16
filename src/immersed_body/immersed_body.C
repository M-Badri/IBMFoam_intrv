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
Description
    class for immersed bodies representation.
SourceFiles
    imm_bodies.C
\*---------------------------------------------------------------------------*/
#include "immersed_body.H"
#include "polyMesh.H"
#include "fvCFD.H"
#include "fvMatrices.H"
#include "geometricOneField.H"

#include "interpolationCellPoint.H"
#include "interpolationCell.H"
#include "meshSearch.H"
#include "List.H"
#include "ListOps.H"

#include "OFstream.H"

#include "cellSet.H"
#include "addToRunTimeSelectionTable.H"
#include "PstreamReduceOps.H"

#include "fvcSmooth.H"
#include "fvMeshSubset.H"
#include "solver_info.H" 

#define ORDER 2

using namespace Foam;

//---------------------------------------------------------------------------//
immersed_body::immersed_body
(
    word body_name,
    const Foam::fvMesh& mesh,
    dictionary& IBMFoamDict,
    dictionary& transportProperties,
    label body_id,
    label recomputeM0,
    std::shared_ptr<geom_model> b_geomModel,
    autoPtr<ib_interpol>& ibIntp,
    List<labelList>& cellPoints
)
:
body_name_(body_name),
isActive_(true),
immersedDict_(IBMFoamDict.subDict(body_name_)),
mesh_(mesh),
trans_properties_(transportProperties),
geom_model_(std::move(b_geomModel)),
cellPoints_(cellPoints),
Axis_(vector::one),
AxisOld_(vector::one),
omega_(0.0),
omegaOld_(0.0),
Vel_(vector::zero),
VelOld_(vector::zero),
a_(vector::zero),
alpha_(vector::zero),
totalAngle_(vector::zero),
CoNum_(0.0),
rhoF_(trans_properties_.lookup("rho")),
body_id_(body_id),
update_torq_(false),
bodyOperation_(0),
octreeField_(mesh_.nCells(), 0),
cellToStartInCreateIB_(0),
startSynced_(false),
totRotMatrix_(tensor::I),
sdBasedLambda_(false),
intSpan_(2.0),
charCellSize_(1e3),
refineBuffers_(0),
recomputeM0_(recomputeM0),
t_to_set_static_(-1),
staticContactPost_(vector::zero)
{
    #include "initIB.H"

    InfoH << iB_Info << "Finished body initialization" << endl;
    InfoH << basic_Info << "New bodyID: " << body_id_ << " name: "
        << body_name_ << " rhoS: " << geom_model_->getRhoS()
        << " dC: " << getDC() << endl;
}
//---------------------------------------------------------------------------//
immersed_body::~immersed_body()
{
}
//---------------------------------------------------------------------------//
void immersed_body::create_immersed_body
(
    volScalarField& body,
    volScalarField& refine_F,
    bool synchCreation
)
{
    geom_model_->create_immersed_body(
        body,
        octreeField_,
        cellPoints_
    );

    if(synchCreation)
    {
        syncCreateImmersedBody(body, refine_F);
    }
}
//---------------------------------------------------------------------------//
void immersed_body::syncCreateImmersedBody
(
    volScalarField& body,
    volScalarField& refine_F
)
{
    geom_model_->setOwner();
    InfoH << iB_Info << "body: " << body_id_
        << " owner: " << geom_model_->getOwner() << endl;

    InfoH << iB_Info << "Computing geometrical properties" << endl;
    geom_model_->calculateGeometricalProperties(body);

    computeBodyCoNumber();

    InfoH << iB_Info << "-- body: " << body_id_
        << " current center of mass position: " << geom_model_->getCoM() << endl;

    const List<DynamicLabelList>& surfCells = geom_model_->getSurfaceCellList();
    DynamicLabelList zeroList(surfCells[Pstream::myProcNo()].size(), 0);

    constructRefineField(
        body,
        refine_F,
        surfCells[Pstream::myProcNo()],
        zeroList
    );

    scalarList charCellSizeL(Pstream::nProcs(),1e4);
    forAll (surfCells[Pstream::myProcNo()],sCellI)
    {
        label cellI = surfCells[Pstream::myProcNo()][sCellI];
        charCellSizeL[Pstream::myProcNo()] =
            min(charCellSizeL[Pstream::myProcNo()],
                Foam::pow(mesh_.V()[cellI],0.3333)
            );
    }
    forAll(charCellSizeL,indl)
    {
        if(charCellSizeL[indl] > 5e3)
        {
            charCellSizeL[indl] = -1.0;
        }
    }

    charCellSize_ = gMax(charCellSizeL);
    InfoH << iB_Info << "Body characteristic cell size: "
        << charCellSize_ << endl;
}
//---------------------------------------------------------------------------//
void immersed_body::syncImmersedBodyParralell1
(
    volScalarField& body,
    volScalarField& refine_F
)
{
    geom_model_->setOwner();
    InfoH << iB_Info << "body: " << body_id_
        << " owner: " << geom_model_->getOwner() << endl;

    InfoH << iB_Info << "Computing geometrical properties" << endl;
    geom_model_->calculateGeometricalPropertiesParallel(body);
}
//---------------------------------------------------------------------------//
void immersed_body::syncImmersedBodyParralell2
(
    volScalarField& body,
    volScalarField& refine_F
)
{

    InfoH << iB_Info << "-- body: " << body_id_
        << " current center of mass position: " << geom_model_->getCoM() << endl;

    const List<DynamicLabelList>& surfCells = geom_model_->getSurfaceCellList();
    DynamicLabelList zeroList(surfCells[Pstream::myProcNo()].size(), 0);

    constructRefineField(
        body,
        refine_F,
        surfCells[Pstream::myProcNo()],
        zeroList
    );

    scalarList charCellSizeL(Pstream::nProcs(),1e4);
    forAll (surfCells[Pstream::myProcNo()],sCellI)
    {
        label cellI = surfCells[Pstream::myProcNo()][sCellI];
        charCellSizeL[Pstream::myProcNo()] =
            min(charCellSizeL[Pstream::myProcNo()],
                Foam::pow(mesh_.V()[cellI],0.3333)
            );
    }
    forAll(charCellSizeL,indl)
    {
        if(charCellSizeL[indl] > 5e3)
        {
            charCellSizeL[indl] = -1.0;
        }
    }

    charCellSize_ = gMax(charCellSizeL);
    InfoH << iB_Info << "Body characteristic cell size: "
        << charCellSize_ << endl;
}
//---------------------------------------------------------------------------//
void immersed_body::constructRefineField
(
    volScalarField& body,
    volScalarField& refine_F,
    DynamicLabelList cellsToIterate,
    DynamicLabelList startLevel
)
{
    if(refineBuffers_ == 0)
        return;

    DynamicLabelList cellsToIterateC;
    DynamicLabelList cellsToIterateF;

    List<DynamicLabelList> cellsToSendToProcs;
    cellsToSendToProcs.setSize(Pstream::nProcs());
    List<DynamicLabelList> cellsToSendToProcsLevel;
    cellsToSendToProcsLevel.setSize(Pstream::nProcs());

    for(label i = 0; i < refineBuffers_; i++)
    {
        forAll(cellsToIterate, cellI)
        {
            if(startLevel[cellI] == i)
                cellsToIterateC.append(cellsToIterate[cellI]);
        }

        forAll(cellsToIterateC, cellI)
        {
            labelList cellFaces(mesh_.cells()[cellsToIterateC[cellI]]);
            forAll(cellFaces, faceI)
            {
                if (mesh_.isInternalFace(cellFaces[faceI]))
                {
                    label nCell(mesh_.owner()[cellFaces[faceI]]);
                    if(nCell == cellsToIterateC[cellI])
                    {
                        nCell = mesh_.neighbour()[cellFaces[faceI]];
                    }

                    if(refine_F[nCell] == 0)
                    {
                        if(i > 0)
                        {
                            if(body[nCell] < SMALL)
                            {
                                refine_F[nCell] = 1;
                                cellsToIterateF.append(nCell);
                            }
                        }
                        else
                        {
                            refine_F[nCell] = 1;
                            cellsToIterateF.append(nCell);
                        }
                    }
                }
                else
                {
                    label facePatchId(mesh_.boundaryMesh().whichPatch(
                        cellFaces[faceI]
                    ));

                    const polyPatch& cPatch = mesh_.boundaryMesh()[facePatchId];
                    if (cPatch.type() == "processor")
                    {
                        const processorPolyPatch& procPatch =
                            refCast<const processorPolyPatch>(cPatch);
                        if (procPatch.myProcNo() == Pstream::myProcNo())
                        {
                            cellsToSendToProcs[procPatch.neighbProcNo()].append(
                                cPatch.whichFace(cellFaces[faceI])
                            );
                            cellsToSendToProcsLevel[procPatch.neighbProcNo()]
                                .append(i+1);
                        }
                        else
                        {
                            cellsToSendToProcs[procPatch.myProcNo()].append(
                                cPatch.whichFace(cellFaces[faceI])
                            );
                            cellsToSendToProcsLevel[procPatch.myProcNo()]
                                .append(i+1);
                        }
                    }
                }
            }
        }
        cellsToIterateC = cellsToIterateF;
        cellsToIterateF.clear();
    }

    List<DynamicLabelList> facesReceivedFromProcs;
    List<DynamicLabelList> cellsReceivedFromProcsLevel;

    // send points that are not on this proc to other proc
    PstreamBuffers pBufs(Pstream::commsTypes::nonBlocking);
    for (label proci = 0; proci < Pstream::nProcs(); proci++)
    {
        if (proci != Pstream::myProcNo())
        {
            UOPstream send(proci, pBufs);
            send << cellsToSendToProcs[proci];
        }
    }
    pBufs.finishedSends();
    // recieve points from other procs
    for (label proci = 0; proci < Pstream::nProcs(); proci++)
    {
        if (proci != Pstream::myProcNo())
        {
            UIPstream recv(proci, pBufs);
            DynamicLabelList recList (recv);
            facesReceivedFromProcs.append(recList);
        }
        else
        {
            DynamicLabelList recList;
            facesReceivedFromProcs.append(recList);
        }
    }
    pBufs.clear();

    for (label proci = 0; proci < Pstream::nProcs(); proci++)
    {
        if (proci != Pstream::myProcNo())
        {

            UOPstream send(proci, pBufs);
            send << cellsToSendToProcsLevel[proci];
        }
    }
    pBufs.finishedSends();
    for (label proci = 0; proci < Pstream::nProcs(); proci++)
    {
        if (proci != Pstream::myProcNo())
        {
            UIPstream recv(proci, pBufs);
            DynamicLabelList recList (recv);
            cellsReceivedFromProcsLevel.append(recList);
        }
        else
        {
            DynamicLabelList recList;
            cellsReceivedFromProcsLevel.append(recList);
        }
    }
    pBufs.clear();


    DynamicLabelList newCellsToIterate;
    DynamicLabelList newCellsToIterateStartLevel;

    for (label otherProci = 0;
        otherProci < facesReceivedFromProcs.size();
        otherProci++
    )
    {
        for (label faceI = 0;
            faceI < facesReceivedFromProcs[otherProci].size();
            faceI++
        )
        {
            label cellProcI(0);
            forAll (mesh_.boundaryMesh(), patchi)
            {
                const polyPatch& cPatch = mesh_.boundaryMesh()[patchi];
                if (cPatch.type() == "processor")
                {
                    const processorPolyPatch& procPatch =
                        refCast<const processorPolyPatch>(cPatch);

                    if (procPatch.myProcNo() == Pstream::myProcNo()
                        && procPatch.neighbProcNo() == otherProci)
                    {
                        cellProcI = mesh_.faceOwner()[cPatch.start()
                            + facesReceivedFromProcs[otherProci][faceI]];
                        if(refine_F[cellProcI] == 0)
                        {
                            newCellsToIterate.append(cellProcI);
                            newCellsToIterateStartLevel.append(
                                cellsReceivedFromProcsLevel[otherProci][faceI]
                            );
                        }
                        break;
                    }
                    else if (procPatch.myProcNo() == otherProci
                        && procPatch.neighbProcNo() == Pstream::myProcNo())
                    {
                        cellProcI = mesh_.faceOwner()[cPatch.start()
                            + facesReceivedFromProcs[otherProci][faceI]];
                        if(refine_F[cellProcI] == 0)
                        {
                            newCellsToIterate.append(cellProcI);
                            newCellsToIterateStartLevel.append(
                                cellsReceivedFromProcsLevel[otherProci][faceI]
                            );
                        }
                        break;
                    }
                }
            }
        }
    }

    bool contBool(false);
    if(newCellsToIterate.size() > 0)
    {
        contBool = true;
    }

    reduce(contBool, orOp<bool>());

    if(contBool)
    {
        constructRefineField
        (
            body,
            refine_F,
            newCellsToIterate,
            newCellsToIterateStartLevel
        );
    }
}
//---------------------------------------------------------------------------//
void immersed_body::post_pimple_update_immersed_body
(
    volScalarField& body,
    volVectorField& f
)
{
    updateCoupling(body,f);

    Vel_ = VelOld_;
    Axis_ = AxisOld_;
    omega_ = omegaOld_;
}
//---------------------------------------------------------------------------//
void immersed_body::updateCoupling
(
    volScalarField& body,
    volVectorField& f
)
{
    vector FV(vector::zero);
    vector TA(vector::zero);

    List<DynamicLabelList> intLists;
    List<DynamicLabelList> surfLists;
    DynamicVectorList refCoMList;

    geom_model_->getReferencedLists(
        intLists,
        surfLists,
        refCoMList
    );

    forAll (intLists, i)
    {
        DynamicLabelList& intListI = intLists[i];
        forAll (intListI, intCell)
        {
            label cellI = intListI[intCell];

            FV -=  f[cellI]*mesh_.V()[cellI];
            TA -=  ((mesh_.C()[cellI] - refCoMList[i])^f[cellI])
                *mesh_.V()[cellI];
        }
    }

    forAll (surfLists, i)
    {
        DynamicLabelList& surfListI = surfLists[i];
        forAll (surfListI, surfCell)
        {
            label cellI = surfListI[surfCell];

            FV -=  f[cellI]*mesh_.V()[cellI];
            TA -=  ((mesh_.C()[cellI] - refCoMList[i])^f[cellI])
                *mesh_.V()[cellI];
        }
    }

  reduce(FV, sumOp<vector>());
  reduce(TA, sumOp<vector>());
  FV *= rhoF_.value();
  TA *= rhoF_.value();

  FCoupling_ = forces(FV, TA);
}
//---------------------------------------------------------------------------//
void immersed_body::updateMovement
(
    scalar deltaT
)
{
    updateMovementComp(deltaT,Vel_,Axis_,omega_);
}
void immersed_body::updateMovement
(
    vector Vel,
    vector Axis,
    scalar omega
)
{
    scalar deltaT = mesh_.time().deltaT().value();
    updateMovementComp(deltaT,Vel,Axis,omega);
}
void immersed_body::updateMovementComp
(
    scalar deltaT,
    vector Vel,
    vector Axis,
    scalar omega
)
{
    auto updateTranslation = [&]()
    {

        const uniformDimensionedVectorField& g =
            mesh_.lookupObject<uniformDimensionedVectorField>("g");
        
        vector FG(vector::zero);
        if(!solver_info::getOnlyDEM())
            FG = geom_model_->getM0()*(1.0-rhoF_.value()
            /geom_model_->getRhoS().value())*g.value();
        else
            FG = geom_model_->getM0()*g.value();

        vector F(FCoupling_.F);
        F += FContact_.F;
        F += FG;

        if(!case3D)
        {
            F[emptyDim] *= 0;
            FG[emptyDim] *= 0;
        }
        if(geom_model_->getM0() > 0)
        {

            InfoH << iB_Info <<"-- body "<< body_id_ <<" ParticelMass  : " << geom_model_->getM0() << endl;
            InfoH << iB_Info <<"-- body "<< body_id_ <<" Acting Force  : " << F << endl;
            a_  = F/(geom_model_->getM0());
            Vel_ = Vel + deltaT*a_;
            InfoH << iB_Info <<"-- body "<< body_id_ <<" accelaration  : " << a_ << endl;
        }
    };

    auto updateRotation = [&]()
    {
        if(mag(geom_model_->getI()) > 0)
        {
            vector T(FCoupling_.T);
            T += FContact_.T;

            alpha_ = inv(geom_model_->getI()) & T;
            vector Omega(Axis*omega + deltaT*alpha_);
            omega_ = mag(Omega);

            if (omega_ < SMALL)
            {
                Axis_ = vector::one;
                if (!case3D)
                {
                    const vector validDirs = (geometricD + vector::one)/2;
                    Axis_ -= validDirs;
                }
            }
            else
            {
                Axis_ =  Omega/(omega_+SMALL);
                if (!case3D)
                {
                    const vector validDirs = (geometricD + vector::one)/2;
                    Axis_ = cmptMultiply(vector::one-validDirs,Axis_);
                }
            }
            Axis_ /= mag(Axis_);
        }
    };

    auto updateRotationFixedAxis = [&]()
    {
        vector T(FCoupling_.T);
        T += FContact_.T;

        vector Omega(Axis*omega + deltaT * (inv(geom_model_->getI()) & T));

        omega_ = mag(Omega);

        vector newAxis = Omega/(omega_+SMALL);
        if ((newAxis & Axis_) < 0) Axis_ *= (-1.0);;
    };

    if (bodyOperation_ == 0 or bodyOperation_ == 3)
    {
        return;
    }
    else if (bodyOperation_ == 1)
    {
        updateRotation();
        return;
    }
    else if (bodyOperation_ == 2)
    {
        updateTranslation();
        return;
    }
    else if (bodyOperation_ == 4)
    {
        updateRotationFixedAxis();
        return;
    }

    updateTranslation();
    if (update_torq_)
    {
        updateRotation();
    }

    return;
}
//---------------------------------------------------------------------------//
void immersed_body::moveImmersedBody
(
    scalar deltaT
)
{
    if (bodyOperation_ == 0) return;


        if (mag(deltaT + 1.0) < SMALL) deltaT = mesh_.time().deltaT().value();


        scalar angle     = omega_*deltaT;

        vector transIncr = Vel_*deltaT;

        tensor rotMatrix(Foam::cos(angle)*tensor::I);
        rotMatrix += Foam::sin(angle)*tensor(
            0.0,      -Axis_.z(),  Axis_.y(),
            Axis_.z(), 0.0,       -Axis_.x(),
            -Axis_.y(), Axis_.x(),  0.0
        );
        rotMatrix += (1.0-Foam::cos(angle))*(Axis_ * Axis_);

        totRotMatrix_ = rotMatrix & totRotMatrix_;
        vector eulerAngles;
        scalar sy = Foam::sqrt(totRotMatrix_.xx()*totRotMatrix_.xx()
            + totRotMatrix_.yy()*totRotMatrix_.yy());

        if (sy > SMALL)
        {
            eulerAngles.x() =
                Foam::atan2(totRotMatrix_.zy(),totRotMatrix_.zz());
            eulerAngles.y() = Foam::atan2(-totRotMatrix_.zx(),sy);
            eulerAngles.z() =
                Foam::atan2(totRotMatrix_.yx(),totRotMatrix_.xx());
        }
        else
        {
            eulerAngles.x() =
                Foam::atan2(-totRotMatrix_.yz(),totRotMatrix_.yy());
            eulerAngles.y() = Foam::atan2(-totRotMatrix_.zx(),sy);
            eulerAngles.z() = 0.0;
        }

        geom_model_->bodyRotatePoints(angle,Axis_);
        geom_model_->bodyMovePoints(transIncr);


  
}
//---------------------------------------------------------------------------//
void immersed_body::printBodyInfo()
{
    InfoH << iB_Info;
    InfoH << "-- body " << body_id_ << " CoM                  : "
        << geom_model_->getCoM() << endl;
    InfoH << "-- body " << body_id_ << " linear velocity      : "
        << Vel_ << endl;
    InfoH << "-- body " << body_id_ << " angluar velocity     : "
        << omega_ << endl;
    InfoH << "-- body " << body_id_ << " axis of rotation     : "
        << Axis_ << endl;
    InfoH << "-- body " << body_id_ << " total rotation matrix: "
        << totRotMatrix_ << endl;
}
//---------------------------------------------------------------------------//
void immersed_body::update_vector_field
(
    volVectorField& VS,
    word VName,
    volScalarField& body
)
{
    word BC = immersedDict_.subDict(VName).lookup("BC");

    List<DynamicLabelList> intLists;
    List<DynamicLabelList> surfLists;
    DynamicVectorList refCoMList;

    geom_model_->getReferencedLists(
        intLists,
        surfLists,
        refCoMList
    );

    if (BC=="noSlip")
    {
        if ( bodyOperation_==0)
        {
            forAll (intLists, i)
            {
                DynamicLabelList& intListI = intLists[i];
                forAll (intListI, intCell)
                {
                    label cellI = intListI[intCell];
                    VS[cellI]   = Vel_;
                }
            }

            forAll (surfLists, i)
            {
                DynamicLabelList& surfListI = surfLists[i];
                forAll (surfListI, surfCell)
                {
                    label cellI = surfListI[surfCell];
                    VS[cellI]   = Vel_;
                }
            }
        }
        else
        {
            forAll (intLists, i)
            {
                DynamicLabelList& intListI = intLists[i];
                forAll (intListI, intCell)
                {
                    label cellI = intListI[intCell];
                    vector planarVec =  mesh_.C()[cellI] - refCoMList[i]
                                    - Axis_*(
                                    (mesh_.C()[cellI] - refCoMList[i])
                                    &Axis_);

                    vector VSvalue = (-(planarVec^Axis_)*omega_ + Vel_);
                    VS[cellI] = VSvalue;
                }
            }

            forAll (surfLists, i)
            {
                DynamicLabelList& surfListI = surfLists[i];
                forAll (surfListI, surfCell)
                {
                    label cellI = surfListI[surfCell];

                    vector planarVec =  mesh_.C()[cellI] - refCoMList[i]
                                    - Axis_*(
                                    (mesh_.C()[cellI] - refCoMList[i])
                                    &Axis_);

                    vector VSvalue = (-(planarVec^Axis_)*omega_ + Vel_);
                    VS[cellI] = VSvalue;
                }
            }
        }
    }
}
vectorField immersed_body::getUatIbPoints()
{
    const List<point>& ibPoints = intpInfo_->getIbPoints();
    vectorField ibPointsVal(ibPoints.size());
    forAll(ibPoints, pointI)
    {

        vector planarVec =  ibPoints[pointI] - geom_model_->getCoM()
                            - Axis_*(
                            (ibPoints[pointI]-geom_model_->getCoM())&Axis_);

        vector VSvalue = (-(planarVec^Axis_)*omega_ + Vel_);
        ibPointsVal[pointI] = VSvalue;
    }

    return ibPointsVal;
}
//---------------------------------------------------------------------------//

void immersed_body::postContactUpdateBodyField
(
    volScalarField& body,
    volScalarField& refine_F
)
{
    geom_model_->resetBody(body);

    create_immersed_body(body,refine_F,false);
}
//---------------------------------------------------------------------------//
void immersed_body::recreateBodyField
(
    volScalarField& body,
    volScalarField& refine_F
)
{
    octreeField_ = Field<label>(mesh_.nCells(), 0);
    geom_model_->getSurfaceCellList()[Pstream::myProcNo()].clear();
    geom_model_->getInternalCellList()[Pstream::myProcNo()].clear();
    create_immersed_body(body,refine_F,false);
}
//---------------------------------------------------------------------------//
void immersed_body::computeBodyCoNumber()
{
    label auxCntr(0);
    scalar VelMag(mag(Vel_));

    meanCoNum_ = 0.0;

    scalar rotCoNumB(omega_*getDC()*0.5*mesh_.time().deltaT().value());

    const List<DynamicLabelList>& surfCells = geom_model_->getSurfaceCellList();
    forAll (surfCells[Pstream::myProcNo()],sCellI)
    {
        label cellI(surfCells[Pstream::myProcNo()][sCellI]);

        scalar dCell(Foam::pow(mesh_.V()[cellI],0.3333));
        scalar CoNumCell(VelMag*mesh_.time().deltaT().value()/dCell);

        CoNumCell+=rotCoNumB/dCell;

        CoNum_      = max(CoNum_,CoNumCell);
        meanCoNum_ += CoNumCell;
        auxCntr    += 1;
    }
    const List<DynamicLabelList>& intCells = geom_model_->getInternalCellList();
    forAll (intCells[Pstream::myProcNo()],iCellI)
    {
        label cellI(intCells[Pstream::myProcNo()][iCellI]);

        scalar dCell(Foam::pow(mesh_.V()[cellI],0.3333));
        scalar CoNumCell(VelMag*mesh_.time().deltaT().value()/dCell);

        CoNum_      = max(CoNum_,CoNumCell);
        meanCoNum_ += CoNumCell;
        auxCntr    += 1;
    }

    reduce(meanCoNum_, sumOp<scalar>());
    reduce(auxCntr, sumOp<scalar>());
    reduce(CoNum_, maxOp<scalar>());

    if(auxCntr > 0)
    {
        meanCoNum_ /= auxCntr;
    }

    InfoH << iB_Info << "-- body " << body_id_
        << " Courant Number mean: " << meanCoNum_
        << " max: " << CoNum_ << endl;

}

//---------------------------------------------------------------------------//
void immersed_body::printMomentum()
{
    vector L(geom_model_->getI()&(Axis_*omega_));
    vector p(geom_model_->getM()*Vel_);

    InfoH << iB_Info;
    InfoH << "-- body " << body_id_ << "  linear momentum:" << p
         << " magnitude: " << mag(p) <<endl;
    InfoH << "-- body " << body_id_ << " angular momentum:" << L
         << " magnitude: " << mag(L) <<endl;
}
//---------------------------------------------------------------------------//
void immersed_body::printStats()
{
    vector L(geom_model_->getI()&(Axis_*omega_));
    vector p(geom_model_->getM()*Vel_);

    InfoH << iB_Info << "-- body " << body_id_ << "  linear momentum:" << p
        << " magnitude: " << mag(p) <<endl;
    InfoH << "-- body " << body_id_ << " angular momentum:" << L
        << " magnitude: " << mag(L) <<endl;
    InfoH << basic_Info << "-- body " << body_id_ << " CoM :"
        << geom_model_->getCoM() << endl;
    InfoH << basic_Info << "-- body " << body_id_ << "  linear velocity:"
        << Vel_ << " magnitude: " << mag(Vel_) <<endl;
    InfoH << "-- body " << body_id_ << " angular velocity:" << omega_
        << " magnitude: " << mag(omega_) <<endl;
    InfoH << "-- body " << body_id_ << "    rotation axis:" << Axis_
        << " magnitude: " << mag(Axis_) <<endl;
}
//---------------------------------------------------------------------------//
void immersed_body::switchActiveOff
(
    volScalarField& body
)
{
    isActive_ = false;
    geom_model_->resetBody(body);
}
//---------------------------------------------------------------------------//
void immersed_body::initSyncWithFlow(const volVectorField& U)
{
    volVectorField curlU(fvc::curl(U));

    vector meanV(vector::zero);
    scalar totVol(0);
    vector meanC(vector::zero);
    label  cellI;
    const List<DynamicLabelList>& intCells = geom_model_->getInternalCellList();
    forAll (intCells[Pstream::myProcNo()],iCellI)
    {
        cellI   = intCells[Pstream::myProcNo()][iCellI];
        meanV  += U[cellI]*mesh_.V()[cellI];
        meanC  += curlU[cellI]*mesh_.V()[cellI];
        totVol += mesh_.V()[cellI];
    }
    reduce(meanV, sumOp<vector>());
    reduce(meanC, sumOp<vector>());
    reduce(totVol, sumOp<scalar>());
    if(totVol > 0)
    {
        Vel_ = meanV/(totVol);
        meanC/=(totVol);
    }
    vector Omega(0.5*meanC);
    if(update_torq_)
    {
        omega_ = mag(Omega);
        if (omega_ < SMALL)
        {
            Axis_ = vector::one;
            if (mesh_.nGeometricD() < 3)
            {
                const vector validDirs =
                    (mesh_.geometricD() + Vector<label>::one)/2;
                Axis_ -= validDirs;
            }
        }
        else
        {
            Axis_ =  Omega/(omega_+SMALL);
            if (mesh_.nGeometricD() < 3)
            {
                const vector validDirs =
                    (mesh_.geometricD() + Vector<label>::one)/2;
                Axis_ = cmptMultiply(vector::one-validDirs,Axis_);
            }
        }
        Axis_ /= mag(Axis_);
    }
    VelOld_     = Vel_;
    omegaOld_   = omega_;
    AxisOld_    = Axis_;
    InfoH << basic_Info << "-- body " << body_id_
        << "initial movement variables:" << endl;
    printStats();
}
//---------------------------------------------------------------------------//
void immersed_body::pimpleUpdate
(
    volScalarField& body,
    volVectorField& f
)
{
    updateCoupling(body, f);
    updateMovement(VelOld_, AxisOld_, omegaOld_);
}
//---------------------------------------------------------------------------//
void immersed_body::check_if_in_domain(volScalarField& body)
{
    if(geom_model_->getM0() < SMALL)
    {
        switchActiveOff(body);
        geom_model_->resetBody(body);
    }

    InfoH << iB_Info << "-- body " << body_id_ << " current M/M0: "
        << geom_model_->getM()/geom_model_->getM0() << endl;
    if (geom_model_->getM()/(geom_model_->getM0()+SMALL) < 1e-2 && case3D)
    {
        switchActiveOff(body);
        geom_model_->resetBody(body);
    }
    else if (!case3D && geom_model_->getNCells() <= 1 && !geom_model_->isCluster())
    {
        switchActiveOff(body);
        geom_model_->resetBody(body);
        InfoH << iB_Info << "-- body " << body_id_ << " switched off" << endl;
    }
}
//---------------------------------------------------------------------------//
void immersed_body::setRestartSim(vector vel, scalar angVel, vector axisRot, bool setStatic, label timesInContact)
{
    Vel_ = vel;
    omega_ = angVel;
    Axis_ = axisRot;
    ibContactClass_->setTimeStepsInContWStatic(timesInContact);
    InfoH << iB_Info << "-- body " << body_id_
        << " timeStepsInContWStatic_: "
        << ibContactClass_->getTimeStepsInContWStatic() << endl;
    if(setStatic)
    {
        bodyOperation_ = 0;
        omega_ = 0;
        Vel_ *= 0;
        InfoH << basic_Info << "-- body " << body_id_ << " set as Static" << endl;
    }
}
//---------------------------------------------------------------------------//
void immersed_body::chceckBodyOp()
{
    if(bodyOperation_ != 5 || t_to_set_static_ == -1)
        return;

    if(!ibContactClass_->checkInContactWithStatic() && ibContactClass_->getTimeStepsInContWStatic() > 0)
    {
        ibContactClass_->setTimeStepsInContWStatic(0);
        return;
    }

    if(ibContactClass_->checkInContactWithStatic())
    {
        ibContactClass_->setTimeStepsInContWStatic(ibContactClass_->getTimeStepsInContWStatic() + 1);
        InfoH << iB_Info << "-- body " << body_id_
            << " timeStepsInContWStatic_: "
            << ibContactClass_->getTimeStepsInContWStatic() << endl;

        if(ibContactClass_->getTimeStepsInContWStatic() == 1)
        {
            staticContactPost_ = geom_model_->getCoM();
        }
        else
        {
            if(mag(staticContactPost_ - geom_model_->getCoM()) > 0.05 * geom_model_->getDC())
            {
                ibContactClass_->setTimeStepsInContWStatic(0);
                return;
            }
        }

        if(ibContactClass_->getTimeStepsInContWStatic() >= t_to_set_static_)
        {
            bodyOperation_ = 0;
            omega_ = 0;
            Vel_ *= 0;
            InfoH << basic_Info << "-- body " << body_id_ << " set as Static" << endl;
        }
    }

    ibContactClass_->inContact_with_static(false);
}
