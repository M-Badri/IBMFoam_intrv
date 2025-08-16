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
    virtual_mesh.C
\*---------------------------------------------------------------------------*/

#ifndef sub_volumeProperties_H
#define sub_volumeProperties_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
/*---------------------------------------------------------------------------*\
                        Class ibContactClass Declaration
\*---------------------------------------------------------------------------*/
class sub_volumeProperties
{
    public:

    sub_volumeProperties(point center):
    center(center)
    {}

    const point center;
    bool isCBody=false;
    bool isTBody=false;
    bool toCheck=true;
    bool isOnEdge=false;

    // point pointInMesh;
    // autoPtr<point> pointInMesh;

    // List<vector> faceNeighbours;
    // autoPtr<List<vector>> faceNeighbours;

    bool isContact()
    {
        return isCBody&&isTBody;
    };

    bool isInBody()
    {
        return  (isCBody||isTBody);
    }
    // Destructors -------------------------------------------------
    ~sub_volumeProperties(){}
};
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
} // End namespace Foam
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
#endif
// ************************************************************************* //
