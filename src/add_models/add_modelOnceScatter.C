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
#include "add_modelOnceScatter.H"
#include "meshSearch.H"

using namespace Foam;

//---------------------------------------------------------------------------//
add_modelOnceScatter::add_modelOnceScatter
(
    const dictionary& add_modelDict,
    const Foam::fvMesh& mesh,
    const bool startTime0,
    std::unique_ptr<geom_model> b_geomModel,
    List<labelList>& cellPoints
)
:
add_model(mesh, std::move(b_geomModel), cellPoints),
add_modelDict_(add_modelDict),
addMode_(word(add_modelDict_.lookup("add_model"))),
bodyAdded_(false),
finishedAddition_(false),

coeffsDict_(add_modelDict_.subDict(addMode_+"Coeffs")),

addDomain_(word(coeffsDict_.lookup("addDomain"))),
scalingMode_(word(coeffsDict_.lookup("scalingMode"))),
rotationMode_(word(coeffsDict_.lookup("rotationMode"))),
addModeI_(word(coeffsDict_.lookup("addMode"))),

addDomainCoeffs_(coeffsDict_.subDict(addDomain_ + "Coeffs")),
scalingModeCoeffs_(coeffsDict_.subDict(scalingMode_ + "Coeffs")),
rotationModeCoeffs_(coeffsDict_.subDict(rotationMode_ + "Coeffs")),
addModeICoeffs_(coeffsDict_.subDict(addModeI_ + "Coeffs")),

partPerAdd_(0),
fieldValue_(0),
addedOnTimeLevel_(0),
startTimeValue_(mesh_.time().value()),

zoneName_(),
minBound_(vector::zero),
maxBound_(vector::zero),

scaleParticles_(false),
minScale_(0),
maxScale_(0),
minScaleFit_(0),
scaleStep_(0),
nTriesBeforeScaling_(0),

rotateParticles_(false),
randomAxis_(false),
axisOfRot_(vector::zero),

bodyAdditionAttemptCounter_(0),
scaleCorrectionCounter_(0),

scaleApplication_(false),
scaleRandomApplication_(false),
rescaleRequirement_(false),
succesfulladition_(false),
scalingFactor_(0),
restartPartCountTemp_(false),
reapeatedAddition_(false),
firstTimeRunning_(true),
cellZoneActive_(false),
boundBoxActive_(false),
octreeField_(mesh_.nCells(), 0),
multiBody_(false),
fieldBased_(false),
fieldCurrentValue_(0),
allActiveCellsInMesh_(true),
randGen_(clock::getTime())
{
    if(!startTime0)
    {
        finishedAddition_ = true;
    }
	init();
}

add_modelOnceScatter::~add_modelOnceScatter()
{
}

//---------------------------------------------------------------------------//
void add_modelOnceScatter::init()
{
    
    cellsInBoundBox_.setSize(Pstream::nProcs());
    cellZonePoints_.setSize(Pstream::nProcs());


	if (addModeI_ == "multiBody")
	{
		partPerAdd_ = (readLabel(addModeICoeffs_.lookup("partPerAdd")));
        multiBody_  = true;
        InfoH << add_model_Info << "-- add_modelMessage-- "
            << "add_model will initialize simulation with "
            << partPerAdd_ << " randomly scattered bodies" << endl;
	}
	else if (addModeI_ == "fieldBased")
	{
		fieldValue_ = (readScalar(addModeICoeffs_.lookup("fieldValue")));
        fieldBased_ = true;
        InfoH << add_model_Info << "-- add_modelMessage-- "
            << "add_model will fill domain up to preset volume fraction" << endl;
		InfoH << "-- add_modelMessage-- "
            << "preset volume fraction: " << fieldValue_ << endl;
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

    // check, if the whole zone is in the mesh
    scalarList procZoneVols(Pstream::nProcs());
    procZoneVols[Pstream::myProcNo()] = 0;
    forAll (cellsInBoundBox_[Pstream::myProcNo()],cellI)
    {
        procZoneVols[Pstream::myProcNo()]+=mesh_.V()[cellsInBoundBox_[Pstream::myProcNo()][cellI]];
    }
    scalar zoneVol(gSum(procZoneVols));
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

	if (scalingMode_ == "noScaling")
	{
		scaleParticles_ = false;
        InfoH << add_model_Info << "-- add_modelMessage-- "
            << "all particles will have the same scale" << endl;
	}
	else if (scalingMode_ == "randomScaling")
	{
		minScale_               = (readScalar(scalingModeCoeffs_.lookup("minScale")));
		maxScale_               = (readScalar(scalingModeCoeffs_.lookup("maxScale")));
		scaleParticles_         = false;
		scaleRandomApplication_ = true;
        InfoH << add_model_Info << "-- add_modelMessage-- "
            << "particles will be randomly scaled" << endl;
	}
	else if (scalingMode_ == "scaleToFit")
	{
		minScaleFit_        = (readScalar(scalingModeCoeffs_.lookup("minScale")));
		scaleStep_          = (readScalar(scalingModeCoeffs_.lookup("scaleStep")));
		nTriesBeforeScaling_= (readScalar(scalingModeCoeffs_.lookup("nTriesBeforeScaling")));
        InfoH << add_model_Info << "-- add_modelMessage-- "
            << "particles will be downscaled to better fill the domain" << endl;
		InfoH << "-- add_modelMessage-- " << "nTriesBeforeDownScaling: "
            << nTriesBeforeScaling_ << endl;
		scaleParticles_ = true;
	}
    else
    {
		InfoH << add_model_Info << "-- add_modelMessage-- "
            << "notImplemented, will crash" << endl;
	}

	if (rotationMode_ == "noRotation")
	{
		rotateParticles_ = false;
		randomAxis_      = false;
        InfoH << add_model_Info << "-- add_modelMessage-- "
            << "source STL will not be rotated upon addition" << endl;
	}
	else if (rotationMode_ == "randomRotation")
	{
		rotateParticles_ = true;
		randomAxis_      = true;
        InfoH << add_model_Info << "-- add_modelMessage-- "
            << "source STL will be randomly rotated upon addition" << endl;
	}
	else if (rotationMode_ == "fixedAxisRandomRotation")
	{
		axisOfRot_       = (rotationModeCoeffs_.lookup("axis"));
        InfoH << add_model_Info << "-- add_modelMessage-- "
            << "source STL will be rotated by a random angle around a fixed axis upon addition" << endl;
		InfoH << "-- add_modelMessage-- " << "set rotation axis: " << axisOfRot_ << endl;
		rotateParticles_ = true;
		randomAxis_      = false;
	}
    else
    {
		InfoH << add_model_Info << "-- add_modelMessage-- "
            << "notImplemented, will crash" << endl;
	}
}

//---------------------------------------------------------------------------//
bool add_modelOnceScatter::shouldAddBody(const volScalarField& body)
{
    if(startTimeValue_ != mesh_.time().value())
    {
        return false;
    }

    if (!finishedAddition_)
    {
        if (multiBody_)
        {

            InfoH << add_model_Info << "-- add_modelMessage-- "
                << "Number of added bodies: " << addedOnTimeLevel_ << endl;

            if (partPerAdd_ < addedOnTimeLevel_)
            {
                finishedAddition_ = false;
            }
            else
            {
                finishedAddition_ = true;
            }
        }

        if (fieldBased_)
        {
            scalar currentLambdaFrac(checkLambdaFraction(body));
            if (currentLambdaFrac < fieldValue_ )
            {
                InfoH << add_model_Info << "-- add_modelMessage-- "
                    << "Current lambda fraction = " << currentLambdaFrac
                    << " < then preset lambda fraction = "
                    << fieldValue_ << endl;
                finishedAddition_ = false;
            }
            else
            {
                finishedAddition_ = true;
            }
        }
    }

    return !finishedAddition_;

}
//---------------------------------------------------------------------------//
std::shared_ptr<geom_model> add_modelOnceScatter::addBody
(
    const   volScalarField& body,
    PtrList<immersed_body>& imm_bodies
)
{
    geom_model_->resetBody();

    bodyAdditionAttemptCounter_++;

    // rotate
    if (rotateParticles_)
    {
        scalar rotAngle = returnRandomAngle();
        if (randomAxis_)
        {
            axisOfRot_ = returnRandomRotationAxis();
        }
        InfoH << add_model_Info << "-- add_modelMessage-- "
            << "Will rotate by " << rotAngle
            << " PiRad around axis " << axisOfRot_ << endl;

        geom_model_->bodyRotatePoints(rotAngle,axisOfRot_);
    }

    // scale
    if (scaleApplication_ or scaleRandomApplication_)
    {
        if (scaleRandomApplication_){scaleStep_ = returnRandomScale();}
        geom_model_->bodyScalePoints(scaleStep_);
    }

    geom_model_->bodyScalePoints(1.02);

    vector CoM(geom_model_->getCoM());
    point bBoxCenter = cellZoneBounds_.midpoint();
    geom_model_->bodyMovePoints(bBoxCenter - CoM);

    vector randomTrans = geom_model_->add_modelReturnRandomPosition(allActiveCellsInMesh_,cellZoneBounds_,randGen_);
    geom_model_->bodyMovePoints(randomTrans);

    // check if the body can be added
    volScalarField helpBodyField_ = body;
    geom_model_->create_immersed_body(
        helpBodyField_,
        octreeField_,
        cellPoints_
    );

    bool canAddBodyI = !isBodyInContact(imm_bodies);

    geom_model_->bodyScalePoints(1.0/1.02);

    reduce(canAddBodyI, andOp<bool>());
    bodyAdded_ = (canAddBodyI);

    if(!bodyAdded_)
	{
		scaleCorrectionCounter_++;
	}

	if(bodyAdded_)
	{
		if(multiBody_)
		{
			InfoH << add_model_Info << "-- add_modelMessage-- "
                << "addedOnTimeLevel:  " << addedOnTimeLevel_<< endl;
			addedOnTimeLevel_++;
		}

		scaleCorrectionCounter_ = 0;

	}

	InfoH << add_model_Info << "-- add_modelMessage-- "
        << "bodyAdditionAttemptNr  : " << bodyAdditionAttemptCounter_<< endl;
	InfoH << "-- add_modelMessage-- "
        << "sameScaleAttempts      : " << scaleCorrectionCounter_<< endl;

	if(scaleCorrectionCounter_ > nTriesBeforeScaling_ && scaleParticles_)
	{
		scaleApplication_ = true;
		rescaleRequirement_ = true;
		scalingFactor_++;
		scaleStep_ = pow(scaleStep_,scalingFactor_);
		if(scaleStep_<minScaleFit_)
		{
			scaleStep_=minScaleFit_;
		}
		scaleCorrectionCounter_ = 0;
	}

    return geom_model_->getCopy();
}
// MODEL SPECIFIC FUNCTIONS==================================================//
//---------------------------------------------------------------------------//
void add_modelOnceScatter::initializeCellZone()
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
void add_modelOnceScatter::updateCellZoneBoundBox()
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
void add_modelOnceScatter::initializeBoundBox()
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
labelList add_modelOnceScatter::getBBoxCellsByOctTree
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
scalar add_modelOnceScatter::checkLambdaFraction(const volScalarField& body)
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
scalar add_modelOnceScatter::returnRandomAngle()
{
    scalar ranNum = 2.0*randGen_.scalar01() - 1.0;
    scalar angle  = ranNum*Foam::constant::mathematical::pi;
	return angle;
}
//---------------------------------------------------------------------------//
scalar add_modelOnceScatter::returnRandomScale()
{
	scalar ranNum       = randGen_.scalar01();
	scalar scaleDiff    = maxScale_ - minScale_;
    scalar scaleFactor  = minScale_ + ranNum*scaleDiff;
	InfoH << add_model_Info << "-- add_modelMessage-- "
        <<"random scaleFactor " << scaleFactor <<endl;
	return scaleFactor;
}
//---------------------------------------------------------------------------//
vector add_modelOnceScatter::returnRandomRotationAxis()
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
