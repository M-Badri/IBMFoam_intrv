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
    class for storage of particle-wall contact variables
\*---------------------------------------------------------------------------*/

#ifndef wallContactVars_H
#define wallContactVars_H

#include "physicalProperties.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{


class wallContactVars
{
    public:
    vector contactCenter_;
    scalar contactVolume_;
    vector contactNormal_;
    scalar contactArea_;

    physicalProperties physicalProperties_;

    vector lVec_;
    vector Veli_;
    scalar Vn_;
    scalar Lc_;
    scalar curAdhN_;    
    
    vector lVector_;    

    vector FtPrev_;     

    // Constructors ------------------------------------------------
    wallContactVars()
    {
        contactCenter_ = vector::zero;
        contactVolume_ = 0;
        contactNormal_ = vector::zero;
        contactArea_ = 0;

        lVec_ = vector::zero;
        Veli_ = vector::zero;
        Vn_ = 0;
        Lc_ = 0;
        curAdhN_ = 0;
        FtPrev_ = vector::zero;
    }

    // Destructors -------------------------------------------------
    ~wallContactVars(){}

    void setMeanCntPars
    (
        const fvMesh&   mesh,
        DynamicList<Tuple2<label,string>>& contactFaces,
        HashTable<physicalProperties,string,Hash<string>>& wallMeanPars
    );

    void setMeanCntPars_Plane
    (
        List<scalar>& contactAreas,
        List<string> contactFaces,
        HashTable<physicalProperties,string,Hash<string>>& wallMeanPars
    );

    physicalProperties& getMeanCntPar()
    {
        return physicalProperties_;
    }

    void clear()
    {
        contactCenter_ = vector::zero;
        contactVolume_ = 0;
        contactNormal_ = vector::zero;
        contactArea_ = 0;

        lVec_ = vector::zero;
        Veli_ = vector::zero;
        Vn_ = 0;
        Lc_ = 0;
        curAdhN_ = 0;
        FtPrev_ = vector::zero;
    }    
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
