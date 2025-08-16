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
#include "add_modelOnceFromFile.H"

#include <memory>

using namespace Foam;

//---------------------------------------------------------------------------//
add_modelOnceFromFile::add_modelOnceFromFile
(
    const dictionary& add_modelDict,
    const Foam::fvMesh& mesh,
    const bool startTime0,
    std::unique_ptr<geom_model> b_geomModel,
    List<labelList>& cellPoints,
    word& b_geom,
    scalar thrSurf
)
:
add_model(mesh, std::move(b_geomModel), cellPoints),
add_modelDict_(add_modelDict),
addMode_(word(add_modelDict_.lookup("add_model"))),
coeffsDict_(add_modelDict_.subDict(addMode_+"Coeffs")),
bodyAdded_(false),
fileName_("constant/" + (word(coeffsDict_.lookup("fileName")))),
ifStream_(fileName_.toAbsolute()),
b_geom_(b_geom),
thrSurf_(thrSurf)
{
    if(!ifStream_.opened())
    {
        FatalErrorIn("add_modelOnceFromFile::add_modelOnceFromFile()")
        << "Can not open IFstream for file: "
        << fileName_.toAbsolute() << nl << exit(FatalError);
    }

    if(!startTime0)
    {
        ifStream_.setEof();
    }
}

add_modelOnceFromFile::~add_modelOnceFromFile()
{
}
//---------------------------------------------------------------------------//
void add_modelOnceFromFile::addSphere(string& line)
{
    IStringStream stringStream(line);
    vector position(stringStream);

    geom_model_->bodyMovePoints(position - geom_model_->getCoM());
}
//---------------------------------------------------------------------------//
void add_modelOnceFromFile::addSTL(string& line)
{
    if(b_geom_ == "convex")
    {
        word stlPath("constant/triSurface/" + line);
        geom_model_ = std::unique_ptr<convex_body>
            (new convex_body(mesh_, stlPath, thrSurf_));
    }
    else
    {
        word stlPath("constant/triSurface/" + line);
        geom_model_ = std::unique_ptr<nonConvex_body>
            (new nonConvex_body(mesh_, stlPath, thrSurf_));
    }
}
//---------------------------------------------------------------------------//
std::shared_ptr<geom_model> add_modelOnceFromFile::addBody
(
    const volScalarField& body,
    PtrList<immersed_body>& imm_bodies
)
{
    string line;
    ifStream_.getLine(line);
    if(line == "")
    {
        InfoH << add_model_Info << "-- add_modelMessage-- "
            << "Skipping empty line at " << ifStream_.lineNumber() - 1 <<endl;
        bodyAdded_ = false;
        return geom_model_->getCopy();
    }

    if (b_geom_ == "sphere")
    {
        addSphere(line);
    }
    else
    {
        addSTL(line);
    }

    bodyAdded_ = true;
    return geom_model_->getCopy();
}
