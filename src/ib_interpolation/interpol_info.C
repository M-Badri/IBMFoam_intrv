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
#include "interpol_info.H"

using namespace Foam;

//---------------------------------------------------------------------------//
interpol_info::interpol_info
(
    const Foam::fvMesh& mesh,
    std::shared_ptr<geom_model>& gModel
)
:
mesh_(mesh),
geom_model_(gModel)
{}

interpol_info::~interpol_info()
{}
//---------------------------------------------------------------------------//