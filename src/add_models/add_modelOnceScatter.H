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
    class for addition and removal of imm_bodies

SourceFiles
    add_modelOnceScatter.C

\*---------------------------------------------------------------------------*/

#ifndef add_modelOnceScatter_H
#define add_modelOnceScatter_H

#include "add_model.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class add_modelOnceScatter
:
    public add_model
{

    private:
    typedef DynamicList<label>  DynamicLabelList;
    const dictionary&           add_modelDict_;
    word                        addMode_;
    bool                        bodyAdded_;
    bool                        finishedAddition_;

    const dictionary&   coeffsDict_;

	word			    addDomain_;
    word			    scalingMode_;
    word			    rotationMode_;
	word			    addModeI_;

	const dictionary& addDomainCoeffs_;
	const dictionary& scalingModeCoeffs_;
	const dictionary& rotationModeCoeffs_;
	const dictionary& addModeICoeffs_;

    label		partPerAdd_;
	scalar		fieldValue_;
    label		addedOnTimeLevel_;
    const scalar startTimeValue_;

	word		zoneName_;
    vector		minBound_;
    vector		maxBound_;

	bool 		scaleParticles_;
	scalar 		minScale_;
	scalar		maxScale_;
	scalar 		minScaleFit_;
	scalar		scaleStep_;
	scalar		nTriesBeforeScaling_;

	bool 		rotateParticles_;
	bool		randomAxis_;
    vector 		axisOfRot_;

	int 		bodyAdditionAttemptCounter_;
	int 		scaleCorrectionCounter_;

	bool 		scaleApplication_;
	bool 		scaleRandomApplication_;
	bool 		rescaleRequirement_;
	bool 		succesfulladition_;
	scalar		scalingFactor_;
	bool		restartPartCountTemp_;
	bool		reapeatedAddition_;
	bool		firstTimeRunning_;
	bool 		cellZoneActive_;
	bool		boundBoxActive_;
	Field<label>octreeField_;
	bool		multiBody_;
	bool		fieldBased_;
	scalar		fieldCurrentValue_;

    List<DynamicLabelList>  cellsInBoundBox_;
    List<pointField>        cellZonePoints_;
    boundBox                cellZoneBounds_;

    bool        allActiveCellsInMesh_;

    Random          randGen_;

    public:

    // Constructors ------------------------------------------------
    add_modelOnceScatter
    (
        const dictionary& add_modelDict,
        const Foam::fvMesh& mesh,
        const bool startTime0,
        std::unique_ptr<geom_model> b_geomModel,
        List<labelList>& cellPoints
    );

    // Destructors -------------------------------------------------
    ~add_modelOnceScatter();

   
    void init();

    bool shouldAddBody(const volScalarField& body);

    std::shared_ptr<geom_model> addBody(
        const volScalarField& body,
        PtrList<immersed_body>& imm_bodies
    );

	scalar returnRandomAngle();

	scalar returnRandomScale();

    vector returnRandomRotationAxis();

	void initializeCellZone();

	void updateCellZoneBoundBox();

    void initializeBoundBox();

	labelList getBBoxCellsByOctTree
    (
        label cellToCheck,
        bool& insideBB,
        vector& bBoxMin,
        vector& bBoxMax,
        List<DynamicLabelList>& bBoxCells
    );

	scalar checkLambdaFraction(const volScalarField& body);

    bool checkIfAdd()
    {
        return true;
    }

    word getAddMode()
    {
        return addMode_;
    }

    bool getBodyAdded()
    {
        return bodyAdded_;
    }

};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
