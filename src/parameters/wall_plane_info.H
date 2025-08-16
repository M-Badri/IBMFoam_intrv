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
    Static class for storing wall material information
\*---------------------------------------------------------------------------*/

#ifndef wall_plane_info_H
#define wall_plane_info_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{


class wall_plane_info
{
    public:

    static const HashTable<List<vector>,string,Hash<string>>& getWallPlaneInfo()
    {
        return wall_plane_info_;
    }

    static void wall_plane_info_insert
    (
        const string& wallPatch,
        const vector& patchNVec,
        const vector& pln_point
    )
    {
        List<vector> planeInfo{{patchNVec,pln_point}};
        wall_plane_info_.insert(wallPatch,planeInfo);
    };

    static HashTable<List<vector>,string,Hash<string>> wall_plane_info_;

};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
