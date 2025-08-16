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
#include "add_modelDistribution.H"
#include "meshSearch.H"

using namespace Foam;

//---------------------------------------------------------------------------//
add_modelDistribution::add_modelDistribution
(
    const dictionary& add_modelDict,
    const Foam::fvMesh& mesh,
    std::unique_ptr<geom_model> b_geomModel,
    List<labelList>& cellPoints
)
:
add_model(mesh, std::move(b_geomModel), cellPoints),
add_modelDict_(add_modelDict),
addMode_(word(add_modelDict_.lookup("add_model"))),
bodyAdded_(false),
distributionDict_
(
    IOobject
    (
        "distributionDict",
        "constant",
        mesh_,
        IOobject::MUST_READ,
        IOobject::NO_WRITE
    )
),
distribution_(scalarList(distributionDict_.lookup("distribution"))),
particleSize_(scalarList(distributionDict_.lookup("particleSize"))),
convertToMeters_(readScalar(distributionDict_.lookup("convertToMeters"))),
volumeOfAddedBodies_(0),

coeffsDict_(add_modelDict_.subDict(addMode_+"Coeffs")),
stlBaseSize_(readScalar(coeffsDict_.lookup("stlBaseSize"))),

addDomain_(word(coeffsDict_.lookup("addDomain"))),
addModeI_(word(coeffsDict_.lookup("addMode"))),

addDomainCoeffs_(coeffsDict_.subDict(addDomain_ + "Coeffs")),
addModeICoeffs_(coeffsDict_.subDict(addModeI_ + "Coeffs")),

useNTimes_(0),
timeBetweenUsage_(0),
partPerAdd_(0),
fieldValue_(0),
addedOnTimeLevel_(0),
partPerAddTemp_(0),

zoneName_(),
minBound_(vector::zero),
maxBound_(vector::zero),

bodyAdditionAttemptCounter_(0),

succesfulladition_(false),
restartPartCountTemp_(false),
reapeatedAddition_(false),
firstTimeRunning_(true),
cellZoneActive_(false),
boundBoxActive_(false),
octreeField_(mesh_.nCells(), 0),
timeBased_(false),
fieldBased_(false),
fieldCurrentValue_(0),
allActiveCellsInMesh_(true),
randGen_(clock::getTime())
{
	init();
}

add_modelDistribution::~add_modelDistribution()
{
}

//---------------------------------------------------------------------------//
void add_modelDistribution::init()
{
    cellsInBoundBox_.setSize(Pstream::nProcs());
    cellZonePoints_.setSize(Pstream::nProcs());
    addedParticlesSize_.setSize(particleSize_.size(), 0);


	if (addModeI_ == "timeBased")
	{
        InfoH << add_model_Info << "-- add_modelMessage-- "
            << "notImplemented, will crash" << endl;
	}
	else if (addModeI_ == "fieldBased")
	{
		fieldValue_ = (readScalar(addModeICoeffs_.lookup("fieldValue")));
        fieldBased_ = true;
        InfoH << add_model_Info << "-- add_modelMessage-- "
            << "add_model will control particles volume fraction" << endl;
		InfoH << "-- add_modelMessage-- " << "preset volume fraction: "
            << fieldValue_ << endl;
	}
    else
    {
        InfoH << add_model_Info << "-- add_modelMessage-- "
            << "notImplemented, will crash" << endl;
    }

	if (addDomain_ == "cellZone")
	{
		zoneName_ = (word(addDomainCoeffs_.lookup("zoneName")));
		cellZoneActive_ = true;
        initializeCellZone();
        InfoH << add_model_Info << "-- add_modelMessage-- "
            << "cellZone based addition zone" << endl;
	}
	else if (addDomain_ == "boundBox")
	{
		minBound_       = (addDomainCoeffs_.lookup("minBound"));
		maxBound_       = (addDomainCoeffs_.lookup("maxBound"));
		boundBoxActive_ = true;

        initializeBoundBox();
        InfoH << add_model_Info << "-- add_modelMessage-- "
            << "boundBox based addition zone" << endl;
	}
	else if (addDomain_ == "domain")
	{
		InfoH << add_model_Info << "-- add_modelMessage-- "
            << "notImplemented, will crash" << endl;
	}
    else
    {
		InfoH << add_model_Info << "-- add_modelMessage-- "
            << "notImplemented, will crash" << endl;
	}

    
    scalarList procZoneVols(Pstream::nProcs());
    procZoneVols[Pstream::myProcNo()] = 0;
    forAll (cellsInBoundBox_[Pstream::myProcNo()],cellI)
    {
        procZoneVols[Pstream::myProcNo()]+=mesh_.V()[cellsInBoundBox_[Pstream::myProcNo()][cellI]];
    }

    Pstream::gatherList(procZoneVols, 0);
    Pstream::scatter(procZoneVols, 0);

    scalar zoneVol(0);
    forAll (procZoneVols, procI)
    {
        zoneVol += procZoneVols[procI];
    }

    scalar zoneBBoxVol(cellZoneBounds_.volume());
    if (zoneVol - zoneBBoxVol > 1e-5*zoneBBoxVol)
    {
        allActiveCellsInMesh_ = false;
        InfoH << add_model_Info << "-- add_modelMessage-- "
             << "addition zone NOT completely immersed in mesh "
             << "this computation will be EXPENSIVE" << endl;
        InfoH << zoneVol << " " << zoneBBoxVol << endl;
    }
    else
    {
        InfoH << add_model_Info << "-- add_modelMessage-- "
             << "addition zone completely immersed in mesh -> OK" << endl;
    }

	partPerAddTemp_ = partPerAdd_;
}

//---------------------------------------------------------------------------//
bool add_modelDistribution::shouldAddBody(const volScalarField& body)
{

    if (timeBased_)
    {
        scalar timeVal(mesh_.time().value());
        scalar deltaTime(mesh_.time().deltaT().value());
        scalar tmFrac(timeVal/timeBetweenUsage_);
        tmFrac -=  floor(tmFrac+deltaTime);

        InfoH << add_model_Info << "-- add_modelMessage-- "
            << "Time/(Time beween usage) - floor(Time/Time beween usage): "
            << tmFrac << endl;

        InfoH << "-- add_modelMessage-- "
            << "Number of bodies added on this time level: "
            << addedOnTimeLevel_ << endl;

        bool tmLevelOk(tmFrac < deltaTime);

        if (not tmLevelOk)
        {
            addedOnTimeLevel_ = 0;
            return false;
        }

        if (partPerAdd_ <= addedOnTimeLevel_) {return false;}

        return (tmLevelOk and useNTimes_ > 0);
    }

    if (fieldBased_)
    {
        scalar currentLambdaFrac(checkLambdaFraction(body));
        if (currentLambdaFrac < fieldValue_ )
        {
            InfoH << add_model_Info << "-- add_modelMessage-- "
                << "Current lambda fraction = " << currentLambdaFrac
                << " < then preset lambda fraction = " << fieldValue_ << endl;
            return true;
        }
    }

    return false;

}
//---------------------------------------------------------------------------//
std::shared_ptr<geom_model> add_modelDistribution::addBody
(
    const volScalarField& body,
    PtrList<immersed_body>& imm_bodies
)
{
    bodyAdditionAttemptCounter_++;

    geom_model_->resetBody();

    Tuple2<label, scalar> scaleFactor = returnScaleFactor();
    InfoH << add_model_Info << "-- add_modelMessage-- "
        << "scaled STL size: " << stlBaseSize_ * scaleFactor.second() << endl;
    geom_model_->bodyScalePoints(scaleFactor.second());
    scalar partVolume(1.0/6.0*3.14*pow(stlBaseSize_ * scaleFactor.second(),3));

    scalar rotAngle = returnRandomAngle();

    vector axisOfRot = returnRandomRotationAxis();

    geom_model_->bodyRotatePoints(rotAngle,axisOfRot);

    vector CoM(geom_model_->getCoM());
    point bBoxCenter = cellZoneBounds_.midpoint();
    geom_model_->bodyMovePoints(bBoxCenter - CoM);

    vector randomTrans = geom_model_->add_modelReturnRandomPosition(allActiveCellsInMesh_,cellZoneBounds_,randGen_);
    geom_model_->bodyMovePoints(randomTrans);
    
    volScalarField helpBodyField_ = body;
    geom_model_->create_immersed_body(
        helpBodyField_,
        octreeField_,
        cellPoints_
    );

    bool canAddBodyI = !isBodyInContact(imm_bodies);

    reduce(canAddBodyI, andOp<bool>());
    bodyAdded_ = (canAddBodyI);

	if(bodyAdded_)
	{
		if(timeBased_)
		{
			InfoH << add_model_Info << "-- add_modelMessage-- "
                << "addedOnTimeLevel:  " << addedOnTimeLevel_<< endl;
			addedOnTimeLevel_++;
			InfoH << add_model_Info << "-- add_modelMessage-- " << "bodyAdded: "
                << bodyAdded_ << " addedOnTimeLevel:  " << addedOnTimeLevel_
                << " useNTimes: " << useNTimes_<<  endl;
			if(addedOnTimeLevel_ == partPerAdd_)
			{
				useNTimes_--;
				InfoH << add_model_Info << "-- add_modelMessage-- "
                    << " useNTimes: " << useNTimes_<<  endl;
				reapeatedAddition_ = false;
			}
		}

		volumeOfAddedBodies_ += partVolume;
        addedParticlesSize_[scaleFactor.first()] += partVolume;
	}

	InfoH << add_model_Info << "-- add_modelMessage-- "
        << "bodyAdditionAttemptNr  : " << bodyAdditionAttemptCounter_<< endl;

    return geom_model_->getCopy();;
}
// MODEL SPECIFIC FUNCTIONS==================================================//
//---------------------------------------------------------------------------//
void add_modelDistribution::initializeCellZone()
{

	label zoneID = mesh_.cellZones().findZoneID(zoneName_);
	InfoH << add_model_Info << "-- add_modelMessage-- "
        << "label of the cellZone " << zoneID << endl;

	const labelList& cellZoneCells = mesh_.cellZones()[zoneID];
    cellsInBoundBox_[Pstream::myProcNo()] = cellZoneCells;

	const pointField& cp = mesh_.C();
	const pointField fCp(cp,cellsInBoundBox_[Pstream::myProcNo()]);
	cellZonePoints_[Pstream::myProcNo()] = fCp;

	updateCellZoneBoundBox();
}
//---------------------------------------------------------------------------//
void add_modelDistribution::updateCellZoneBoundBox()
{
		boundBox cellZoneBounds(cellZonePoints_[Pstream::myProcNo()]);

        reduce(cellZoneBounds.min(), minOp<vector>());
        reduce(cellZoneBounds.max(), maxOp<vector>());

        if (Pstream::myProcNo() == 0)
        {
            minBound_ = cellZoneBounds_.min();
            maxBound_ = cellZoneBounds_.max();
            cellZoneBounds_ = boundBox(minBound_,maxBound_);
        }
}
//---------------------------------------------------------------------------//
void add_modelDistribution::initializeBoundBox()
{
    octreeField_ *= 0;
    List<DynamicLabelList> bBoxCells(Pstream::nProcs());

    bool isInsideBB(false);
    labelList nextToCheck(1,0);
    label iterCount(0);label iterMax(mesh_.nCells());
    while ((nextToCheck.size() > 0 or not isInsideBB) and iterCount < iterMax)
    {
        iterCount++;
        DynamicLabelList auxToCheck;

        forAll (nextToCheck,cellToCheck)
        {
            auxToCheck.append(
                getBBoxCellsByOctTree(
                    nextToCheck[cellToCheck],
                    isInsideBB,
                    minBound_,maxBound_,bBoxCells
                )
            );
        }
        nextToCheck = auxToCheck;
    }

    cellsInBoundBox_[Pstream::myProcNo()] = bBoxCells[Pstream::myProcNo()];

    cellZoneBounds_ = boundBox(minBound_,maxBound_);
}
//---------------------------------------------------------------------------//
void add_modelDistribution::recreateBoundBox()
{
    octreeField_ = Field<label>(mesh_.nCells(), 0);
    cellsInBoundBox_[Pstream::myProcNo()].clear();
    List<DynamicLabelList> bBoxCells(Pstream::nProcs());

    bool isInsideBB(false);
    labelList nextToCheck(1,0);
    label iterCount(0);label iterMax(mesh_.nCells());
    while ((nextToCheck.size() > 0 or not isInsideBB) and iterCount < iterMax)
    {
        iterCount++;
        DynamicLabelList auxToCheck;

        forAll (nextToCheck,cellToCheck)
        {
            auxToCheck.append(
                getBBoxCellsByOctTree(
                    nextToCheck[cellToCheck],
                    isInsideBB,
                    minBound_,maxBound_,bBoxCells
                )
            );
        }
        nextToCheck = auxToCheck;
    }

    cellsInBoundBox_[Pstream::myProcNo()] = bBoxCells[Pstream::myProcNo()];

    cellZoneBounds_ = boundBox(minBound_,maxBound_);
}
//---------------------------------------------------------------------------//
labelList add_modelDistribution::getBBoxCellsByOctTree
(
    label cellToCheck,
    bool& insideBB,
    vector& bBoxMin,
    vector& bBoxMax,
    List<DynamicLabelList>& bBoxCells
)
{
    labelList retList;

    if (octreeField_[cellToCheck] ==0)
    {
        octreeField_[cellToCheck] = 1;
        vector cCenter = mesh_.C()[cellToCheck];
        label   partCheck(0);
        forAll (bBoxMin,vecI)
        {
            if (cCenter[vecI] >= bBoxMin[vecI] and cCenter[vecI] <= bBoxMax[vecI])
            {
                partCheck++;
            }
        }
        bool cellInside = (partCheck == 3) ? true : false;
        if (cellInside)
        {
            bBoxCells[Pstream::myProcNo()].append(cellToCheck);
            insideBB = true;
        }
        if (not insideBB or cellInside)
        {
            retList = mesh_.cellCells()[cellToCheck];
        }
    }
    return retList;
}
//---------------------------------------------------------------------------//
scalar add_modelDistribution::checkLambdaFraction(const volScalarField& body)
{
	scalarList lambdaIntegrate(Pstream::nProcs());
    scalarList volumeIntegrate(Pstream::nProcs());
	scalar lambdaFraction(0);
    forAll (lambdaIntegrate,k)
    {
        lambdaIntegrate[k] = 0;
        volumeIntegrate[k] = 0;
    }
	forAll (cellsInBoundBox_[Pstream::myProcNo()],k)
	{
		label cell = cellsInBoundBox_[Pstream::myProcNo()][k];
		lambdaIntegrate[Pstream::myProcNo()] += mesh_.V()[cell]*body[cell];
		volumeIntegrate[Pstream::myProcNo()] += mesh_.V()[cell];
	}
	lambdaFraction = gSum(lambdaIntegrate)/gSum(volumeIntegrate);
	InfoH << add_model_Info << "-- add_modelMessage-- "
        << "lambda fraction in controlled region: " << lambdaFraction<< endl;
	return lambdaFraction;
}
//---------------------------------------------------------------------------//
scalar add_modelDistribution::returnRandomAngle()
{
    scalar ranNum = 2.0*randGen_.scalar01() - 1.0;
    scalar angle  = ranNum*Foam::constant::mathematical::pi;
	return angle;
}
//---------------------------------------------------------------------------//
vector add_modelDistribution::returnRandomRotationAxis()
{
	vector  axisOfRotation(vector::zero);
	scalar ranNum = 0;

	for (int i=0;i<3;i++)
	{
		ranNum = randGen_.scalar01();
		axisOfRotation[i] = ranNum;
	}

	axisOfRotation /=mag(axisOfRotation);
	return axisOfRotation;
}
//---------------------------------------------------------------------------//
Tuple2<label, scalar> add_modelDistribution::returnScaleFactor()
{
    DynamicLabelList  missingParticles;
    scalar maxSizeVolume = 1.0/6.0*3.14*pow(stlBaseSize_
                           *particleSize_[particleSize_.size()-1]
                           *convertToMeters_/stlBaseSize_,3);
    forAll (addedParticlesSize_,size)
    {
        scalar distribDiff = distribution_[size] - 100*addedParticlesSize_[size]/(volumeOfAddedBodies_+SMALL);
        if(distribDiff > 0)
        {
            scalar meanFactor = particleSize_[size]*convertToMeters_/stlBaseSize_;
            scalar meanVolume = 1.0/6.0*3.14*pow(stlBaseSize_ * meanFactor,3);
            missingParticles.append(floor((distribDiff*maxSizeVolume)/meanVolume));
        }
        else
        {
            missingParticles.append(0);
        }
    }

    label totalMissParts(0);
    forAll (missingParticles,size)
    {
        totalMissParts += missingParticles[size];
    }

    label randomMissPart = floor(randGen_.scalar01()*totalMissParts);
    label missingPart(0);
    forAll (missingParticles,size)
    {
        missingPart = size;
        randomMissPart -= missingParticles[size];
        if(randomMissPart < 0)
        {
            break;
        }
    }

    scalar factor(particleSize_[missingPart - 1] + (particleSize_[missingPart] - particleSize_[missingPart - 1]) * randGen_.scalar01());
    factor *= convertToMeters_/stlBaseSize_;

    Tuple2<label, scalar> returnValue(missingPart, factor);

    return returnValue;
}
