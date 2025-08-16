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
    add_modelDistribution.C
\*---------------------------------------------------------------------------*/

#ifndef add_modelDistribution_H
#define add_modelDistribution_H

#include "add_model.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{


class add_modelDistribution
:
    public add_model
{

    private:
    typedef DynamicList<label>  DynamicLabelList;
    typedef DynamicList<scalar>  DynamicScalarList;
    const dictionary&           add_modelDict_;
    word                        addMode_;
    bool                        bodyAdded_;

    IOdictionary           distributionDict_;
    const scalarList       distribution_;
    const scalarList       particleSize_;
    const scalar       convertToMeters_;

    DynamicScalarList  addedParticlesSize_;
    scalar               volumeOfAddedBodies_;

    const dictionary&   coeffsDict_;
    const scalar       stlBaseSize_;

	word			    addDomain_;
	word			    addModeI_;

	const dictionary& addDomainCoeffs_;
	const dictionary& addModeICoeffs_;

    label       useNTimes_;
    scalar      timeBetweenUsage_;
    label		partPerAdd_;
	scalar		fieldValue_;
    label		addedOnTimeLevel_;
    label		partPerAddTemp_;

	word		zoneName_;
    vector		minBound_;
    vector		maxBound_;

	int 		bodyAdditionAttemptCounter_;
	int 		scaleCorrectionCounter_;

	bool 		succesfulladition_;
	bool		restartPartCountTemp_;
	bool		reapeatedAddition_;
	bool		firstTimeRunning_;
	bool 		cellZoneActive_;
	bool		boundBoxActive_;
	Field<label>octreeField_;
	bool		timeBased_;
	bool		fieldBased_;
	scalar		fieldCurrentValue_;

    List<DynamicLabelList>  cellsInBoundBox_;
    List<pointField>        cellZonePoints_;
    boundBox                cellZoneBounds_;

    bool        allActiveCellsInMesh_;

    Random          randGen_;

    protected:


    public:

    // Constructors ------------------------------------------------
    add_modelDistribution
    (
        const dictionary& add_modelDict,
        const Foam::fvMesh& mesh,
        std::unique_ptr<geom_model> b_geomModel,
        List<labelList>& cellPoints
    );

    // Destructors -------------------------------------------------
    ~add_modelDistribution();

    // Public methods ----------------------------------------------
    void init();

    bool shouldAddBody(const volScalarField& body);

    std::shared_ptr<geom_model> addBody(
        const volScalarField& body,
        PtrList<immersed_body>& imm_bodies
        );

    Tuple2<label, scalar> returnScaleFactor();

	void initializeCellZone();

	void updateCellZoneBoundBox();

    void initializeBoundBox();

    void recreateBoundBox();

	labelList getBBoxCellsByOctTree
    (
        label cellToCheck,
        bool& insideBB,
        vector& bBoxMin,
        vector& bBoxMax,
        List<DynamicLabelList>& bBoxCells
    );

	scalar checkLambdaFraction(const volScalarField& body);

    bool checkIfAdd() {return true;}

    word getAddMode()
    {
        return addMode_;
    }

    bool getBodyAdded()
    {
        return bodyAdded_;
    }

    scalar returnRandomAngle();

    vector returnRandomRotationAxis();
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
