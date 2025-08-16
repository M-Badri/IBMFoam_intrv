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
    add_modelOnceFromFile.C
\*---------------------------------------------------------------------------*/

#ifndef add_modelOnceFromFile_H
#define add_modelOnceFromFile_H

#include "add_model.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class add_modelOnceFromFile
:
    public add_model
{

    private:
    const dictionary&       add_modelDict_;
    word                    addMode_;
    const dictionary&       coeffsDict_;
    bool                    bodyAdded_;

    fileName                fileName_;
    IFstream                ifStream_;

    word                    b_geom_;
    scalar                  thrSurf_;

    // Private Member Functions ------------------------------------

    void addSphere(string& line);

    void addSTL(string& line);

    public:

    // Constructors ------------------------------------------------
    add_modelOnceFromFile
    (
        const dictionary& add_modelDict,
        const Foam::fvMesh& mesh,
        const bool startTime0,
        std::unique_ptr<geom_model> b_geomModel,
        List<labelList>& cellPoints,
        word& b_geom,
        scalar thrSurf
    );

    // Destructors -------------------------------------------------
    ~add_modelOnceFromFile();

    // Public methods ----------------------------------------------

    bool shouldAddBody(const volScalarField& body)
    {
        return not ifStream_.eof();
    }

    std::shared_ptr<geom_model> addBody(
        const volScalarField& body,
        PtrList<immersed_body>& imm_bodies
        );

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
