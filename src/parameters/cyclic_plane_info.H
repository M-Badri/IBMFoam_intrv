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

#ifndef cyclic_plane_info_H
#define cyclic_plane_info_H

#include <map>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{



class cyclic_plane_info
{
    public:

    static const HashTable<List<vector>,string,Hash<string>>& getCyclicPlaneInfo()
    {
        return cyclic_plane_info_;
    }

    static const vector getCyclicTransVec(const string& cyclicPatch)
    {
        const string& cyclicPatchNeighbour = cyclicPatchPairs_[cyclicPatch];

        vector transVec = cyclic_plane_info_[cyclicPatchNeighbour][1] - cyclic_plane_info_[cyclicPatch][1];

        Info << "transVec1: " << transVec << endl;

        transVec = cyclic_plane_info_[cyclicPatch][0] * (transVec & cyclic_plane_info_[cyclicPatch][0]);

        Info << "transVec2: " << transVec << endl;

        return transVec;
    }

    static void insert
    (
        const string& cyclicPatch,
        const vector& patchNVec,
        const vector& pln_point,
        const string& cyclicPatchNeighbour
    )
    {
        List<vector> planeInfo{{patchNVec / mag(patchNVec),pln_point}};
        cyclic_plane_info_.insert(cyclicPatch,planeInfo);

        if (cyclicPatchPairs_.find(cyclicPatch) == cyclicPatchPairs_.end())
        {
            cyclicPatchPairs_.insert(std::pair<string,string>(cyclicPatch,cyclicPatchNeighbour));
        }
    };

    static HashTable<List<vector>,string,Hash<string>> cyclic_plane_info_;
    static std::map<string, string> cyclicPatchPairs_;

};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
