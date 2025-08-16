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
    class for storage of contact variables
\*---------------------------------------------------------------------------*/

#ifndef forces_H
#define forces_H

#include "dictionary.H"
#include "fvCFD.H"

#include "geom_model.H"
#include "ibContactVars.H"
#include "wallContactVars.H"
#include "material_info.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

struct forces{
    vector F = vector::zero;
    vector T = vector::zero;

    forces()
    {}

    forces(const forces& rhs)
    {
        F = rhs.F;
        T = rhs.T;
    }

    forces(const vector& F, const vector& T):
        F(F),
        T(T)
    {}

    void operator=(const forces& rhs)
    {
        F = rhs.F;
        T = rhs.T;
    }

    void operator+=(const forces& rhs)
    {
        F += rhs.F;
        T += rhs.T;
    }

    void operator-=(const forces& rhs)
    {
        F -= rhs.F;
        T -= rhs.T;
    }

    void clear()
    {
        F = vector::zero;
        T = vector::zero;
    }
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
