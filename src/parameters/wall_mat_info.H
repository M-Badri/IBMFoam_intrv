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

#ifndef wall_mat_info_H
#define wall_mat_info_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{


class wall_mat_info
{
    public:

    static const HashTable<material_info,string,Hash<string>>& getWallMatInfo()
    {
        return wall_mat_info_;
    }

    static List<string> getWallPatches()
    {
        return wall_mat_info_.toc();
    }

    static void wall_mat_info_insert
    (
        const string& wallPatch,
        const material_info& matInfo
    )
    {
        wall_mat_info_.insert(wallPatch,matInfo);
    };

    static HashTable<material_info,string,Hash<string>> wall_mat_info_;

};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
