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
#include "add_modelRepeatSamePosition.H"

using namespace Foam;

//---------------------------------------------------------------------------//
add_modelRepeatSamePosition::add_modelRepeatSamePosition
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
coeffsDict_(add_modelDict_.subDict(addMode_+"Coeffs")),
useNTimes_(readLabel(coeffsDict_.lookup("useNTimes"))),
timeBetweenUsage_(readScalar(coeffsDict_.lookup("timeBetweenUsage"))),
addedOnTimeLevel_(0)
{}

add_modelRepeatSamePosition::~add_modelRepeatSamePosition()
{
}


//---------------------------------------------------------------------------//
bool add_modelRepeatSamePosition::shouldAddBody(const volScalarField& body)
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

    if (not tmLevelOk){addedOnTimeLevel_ = 0;}

    return (tmLevelOk and useNTimes_ > 0 and addedOnTimeLevel_ == 0);
}

std::shared_ptr<geom_model> add_modelRepeatSamePosition::addBody
(
    const volScalarField& body,
    PtrList<immersed_body>& imm_bodies
)
{
    volScalarField helpBodyField_ = body;
    geom_model_->create_immersed_body(
        helpBodyField_,
        octreeField_,
        cellPoints_
    );

    bool canAddBodyI = !isBodyInContact(imm_bodies);

    reduce(canAddBodyI, andOp<bool>());

    bodyAdded_ = canAddBodyI;
    if (bodyAdded_) {useNTimes_--;}
    addedOnTimeLevel_++;

    InfoH << add_model_Info << "-- add_modelMessage-- "
        << "will try to use the body " << useNTimes_ << " more times" << endl;

    return geom_model_->getCopy();
}
