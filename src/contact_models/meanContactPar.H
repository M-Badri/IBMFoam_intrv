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
    class for storage of physical properties
\*---------------------------------------------------------------------------*/

#ifndef physicalProperties_H
#define physicalProperties_H

#include "fvCFD.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{


struct physicalProperties
{
    physicalProperties()
    :
        aY_(0),
        aG_(0),
        aMu_(0),
        maxAdhN_(0),
        curAdhN_(0),
        reduceM_(0),
        reduceBeta_(0)
    {}

    physicalProperties
    (
        const scalar aY,
        const scalar aG,
        const scalar aMu,
        const scalar maxAdhN,
        const scalar curAdhN,
        const scalar reduceM,
        const scalar reduceBeta
    ):
        aY_(aY),
        aG_(aG),
        aMu_(aMu),
        maxAdhN_(maxAdhN),
        curAdhN_(curAdhN),
        reduceM_(reduceM),
        reduceBeta_(reduceBeta),
    {}

    // Mean model parame,ters
    scalar aY_ = 0;         
    scalar aG_ = 0;         
    scalar aMu_ = 0;        
    scalar maxAdhN_ = 0;    
    scalar curAdhN_ = 0;    
    scalar reduceM_ = 0;    
    scalar reduceBeta_ = 0;   
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
