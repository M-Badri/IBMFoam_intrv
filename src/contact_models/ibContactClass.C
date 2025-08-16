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
#include "ibContactClass.H"

#include "materialProperties.H"

using namespace Foam;

//---------------------------------------------------------------------------//
ibContactClass::ibContactClass
(
    std::shared_ptr<geom_model>& geom_model,
    const string& material
)
:
geom_model_(geom_model),
isInWallContact_(false),
inContact_with_static_(false),
timeStepsInContWStatic_(0),
matInfo_(materialProperties::getMatProps()[material])
{
}

ibContactClass::ibContactClass(const ibContactClass& other)
:
geom_model_(other.geom_model_),
isInWallContact_(other.isInWallContact_),
inContact_with_static_(other.inContact_with_static_),
timeStepsInContWStatic_(other.timeStepsInContWStatic_),
matInfo_(other.matInfo_)
{
}

ibContactClass::~ibContactClass()
{
}
//---------------------------------------------------------------------------//
