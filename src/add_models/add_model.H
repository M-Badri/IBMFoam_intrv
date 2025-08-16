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
    add_model.C
\*---------------------------------------------------------------------------*/

#ifndef add_model_H
#define add_model_H

#include "dictionary.H"
#include "fvCFD.H"
#include "triSurfaceMesh.H"

#include "geom_model.H"
#include "convex_body.H"
#include "nonConvex_body.H"
#include "sphere_body.H"
#include "periodicBody.H"

#include "contactModel.H"
#include "immersed_body.H"
#include "ibContactVars.H"

#include "virtual_mesh.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class add_model
{

    private:

    protected:

    const  fvMesh&       mesh_;
    List<labelList>& cellPoints_;
    Field<label> octreeField_;
    std::unique_ptr<geom_model> geom_model_;

    public:

    // Constructors ------------------------------------------------
    add_model(
        const fvMesh& mesh,
        std::unique_ptr<geom_model> b_geomModel,
        List<labelList>& cellPoints
    )
    :
    mesh_(mesh),
    cellPoints_(cellPoints),
    geom_model_(std::move(b_geomModel))
    {
        octreeField_ = Field<label>(mesh_.nCells(), 0);
    }

    // Destructors -------------------------------------------------
    virtual ~add_model();

    // -- init
    virtual void init()
    {
        InfoH << add_model_Info << "-- add_modelMessage-- "
            << "notIplemented"<<endl;
    }
    

    virtual bool shouldAddBody(const volScalarField& body)
    {
        InfoH << add_model_Info << "-- add_modelMessage-- "
            << "notIplemented"<<endl;
        return false;
    }
   
    virtual std::shared_ptr<geom_model> addBody(
        const volScalarField& body,
        PtrList<immersed_body>& imm_bodies
        ) = 0;
   

    // -- getters
    virtual word getAddMode()
    {
        InfoH << add_model_Info << "-- add_modelMessage-- "
            << "notIplemented"<<endl;
        return "NO";
    }

    virtual bool getBodyAdded()
    {
        InfoH << add_model_Info << "-- add_modelMessage-- "
            << "notIplemented"<<endl;
        return false;
    }

    virtual void recreateBoundBox()
    {
        InfoH << add_model_Info << "-- add_modelMessage-- "
            << "notIplemented"<<endl;
    }

    bool isBodyInContact(PtrList<immersed_body>& imm_bodies);
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
