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

#ifndef contactVars_H
#define contactVars_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                        Class ibContactVars Declaration
\*---------------------------------------------------------------------------*/

class ibContactVars
{
    public:
    label body_id_;
    vector& Vel_;
    scalar& omega_;
    vector& Axis_;
    scalar& M0_;
    scalar& M_;
    dimensionedScalar rhoS_;

    // Constructors ------------------------------------------------
    ibContactVars
    (
        label body_id,
        vector& Vel,
        scalar& omega,
        vector& Axis,
        scalar& M0,
        scalar& M,
        dimensionedScalar rhoS
    )
    :
    body_id_(body_id),
    Vel_(Vel),
    omega_(omega),
    Axis_(Axis),
    M0_(M0),
    M_(M),
    rhoS_(rhoS)
    {}
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
