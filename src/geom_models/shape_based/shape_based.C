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
#include "shape_based.H"

using namespace Foam;

//---------------------------------------------------------------------------//
shape_based::shape_based
(
    const  fvMesh&   mesh,
    const contactType cType,
    scalar  thrSurf
)
:
geom_model(mesh,cType,thrSurf)
{}
//---------------------------------------------------------------------------//
vector shape_based::add_modelReturnRandomPosition
(
    const bool allActiveCellsInMesh,
    const boundBox  cellZoneBounds,
    Random&          randGen
)
{
    vector ranVec(vector::zero);

    meshSearch searchEng(mesh_);


    vector CoM(getCoM());

    const vector validDirs = (geometricD + vector::one)/2;
    vector dirCorr(cmptMultiply((vector::one - validDirs),CoM));
    dirCorr += cmptMultiply((vector::one - validDirs),0.5*(mesh_.bounds().max() + mesh_.bounds().min()));


    boundBox bodyBounds(getBounds());

    vector maxScales(cellZoneBounds.max() - bodyBounds.max());
    maxScales -= cellZoneBounds.min() - bodyBounds.min();
    maxScales *= 0.5*0.9;

    InfoH << add_model_Info << "-- add_modelMessage-- "
        << "acceptable movements: " << maxScales << endl;

    scalar ranNum = 0;
    for (int i=0;i<3;i++)
    {
        ranNum = 2.0*maxScales[i]*randGen.scalar01() - 1.0*maxScales[i];
        ranVec[i] = ranNum;
    }

    ranVec = cmptMultiply(validDirs,ranVec);
    ranVec += dirCorr;

    return ranVec;
}
//---------------------------------------------------------------------------//
