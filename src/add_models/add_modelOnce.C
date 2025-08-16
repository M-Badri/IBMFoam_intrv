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
#include "add_modelOnce.H"

using namespace Foam;

//---------------------------------------------------------------------------//
add_modelOnce::add_modelOnce
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
bodyAdded_(false)
{
    if(!startTime0)
        bodyAdded_ = true;
}

add_modelOnce::~add_modelOnce()
{
}
//---------------------------------------------------------------------------//

std::shared_ptr<geom_model> add_modelOnce::addBody
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

    bodyAdded_ = true;
    return geom_model_->getCopy();
}
