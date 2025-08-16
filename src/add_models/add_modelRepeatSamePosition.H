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
    add_modelRepeatSamePosition.C

\*---------------------------------------------------------------------------*/

#ifndef add_modelRepeatSamePosition_H
#define add_modelRepeatSamePosition_H

#include "add_model.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{


class add_modelRepeatSamePosition
:
    public add_model
{

    private:
    const dictionary&       add_modelDict_;
    word                    addMode_;
    bool                    bodyAdded_;

    const dictionary&       coeffsDict_;
    label                   useNTimes_;
    scalar                  timeBetweenUsage_;
    label                   addedOnTimeLevel_;

    public:

    // Constructors ------------------------------------------------
    add_modelRepeatSamePosition
    (
        const dictionary& add_modelDict,
        const Foam::fvMesh& mesh,
        std::unique_ptr<geom_model> b_geomModel,
        List<labelList>& cellPoints
    );

    // Destructors -------------------------------------------------
    ~add_modelRepeatSamePosition();

    // Public methods ----------------------------------------------
    //~ void init();

    bool shouldAddBody(const volScalarField& body);

    std::shared_ptr<geom_model> addBody(
        const volScalarField& body,
        PtrList<immersed_body>& imm_bodies
    );

    bool checkIfAdd()
    {
        return true;
    }

    void translateRandom()
    {
        InfoH << add_model_Info << "-- add_modelMessage-- "
            << "Not needed"<< endl;
    }

    void rotateRandom()
    {
        InfoH << add_model_Info << "-- add_modelMessage-- "
            << "Not needed"<< endl;
    }

    word getAddMode()
    {
        return addMode_;
    }

    bool getBodyAdded(){return bodyAdded_;}
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
