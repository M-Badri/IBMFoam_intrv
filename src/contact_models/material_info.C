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
#include "material_info.H"

using namespace Foam;

//---------------------------------------------------------------------------//
material_info::material_info
(
    string material,
    scalar Y,
    scalar nu,
    scalar mu,
    scalar adh_coef,
    scalar eps
)
:
material_(material),
Y_(Y),
nu_(nu),
mu_(mu),
adh_coef_(adh_coef),
eps_(eps)

{
}

material_info::material_info(const material_info& mi)
:
material_(mi.material_),
Y_(mi.Y_),
nu_(mi.nu_),
mu_(mi.mu_),
adh_coef_(mi.adh_coef_),
eps_(mi.eps_)
{
}

material_info::~material_info()
{
}
//---------------------------------------------------------------------------//
