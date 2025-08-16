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
    contactModel

Description
    methods to solve particle-particle contact

SourceFiles
    spectator_mesh.C
\*---------------------------------------------------------------------------*/

#ifndef sMProperties_H
#define sMProperties_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
/*---------------------------------------------------------------------------*\
                        Class ibContactClass Declaration
\*---------------------------------------------------------------------------*/
class sMProperties
{
    public:

    sMProperties(point center):
    center(center)
    {}
    
    const point center;
    point initPoint = vector::zero;
    bool isCBody = false;
    bool isMesh  = false;
    bool toCheck = true;
    bool isInternal = false;
    bool stlFound = false;
    bool initPointSet = false;
    
    // Destructors -------------------------------------------------
    ~sMProperties(){}
};
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
} // End namespace Foam
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
#endif
// ************************************************************************* //
