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

SourceFiles
    wallSubContactInfo.C
\*---------------------------------------------------------------------------*/

#ifndef wallSubContactInfo_H
#define wallSubContactInfo_H

#include "dictionary.H"
#include "fvCFD.H"

#include "geom_model.H"
#include "ibContactVars.H"
#include "ibContactClass.H"
#include "wallContactVars.H"
#include "material_info.H"
#include "physicalProperties.H"
#include "forces.H"
#include "spectator_mesh.H"
#include "virtual_meshWall.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{


class wallSubContactInfo
{
    private:

    typedef DynamicList<label> DynamicLabelList;
    typedef DynamicList<scalar> DynamicScalarList;
    typedef DynamicList<vector> DynamicVectorList;
    typedef DynamicList<point> DynamicPointList;

    List<string> contactPatches_;

   
    List<Tuple2<point,boundBox>> internalBBData_;

   
    HashTable<physicalProperties,string,Hash<string>> wallMeanPars_;

    boundBox BB_;

    label body_id_;

    wallContactVars wallCntvar_;

    scalar reduceM_;    

    forces outForce_; 

    List<autoPtr<virtual_meshWallInfo>> vmWInfoList_; 

    List<autoPtr<virtual_meshWallInfo>> vmPlaneInfoList_;

    bool contact_resolved_ = false;

    public:

    // Constructors ------------------------------------------------
    wallSubContactInfo
    (
    List<Tuple2<point,boundBox>> contactBBData,
    List<Tuple2<point,boundBox>> planeBBData,
    List<string> contactPatches,
    List<Tuple2<point,boundBox>> internalBBData,
    HashTable<physicalProperties,string,Hash<string>> wallMeanPars,
    boundBox BB,
    label body_id
    );

    // Destructors -------------------------------------------------
    ~wallSubContactInfo();

    label getVMContactSize()
    {
        return vmWInfoList_.size();
    }
    
    label getVMPlaneSize()
    {
        return vmPlaneInfoList_.size();
    }

    autoPtr<virtual_meshWallInfo>& getVMContactInfo
    (
        label ID
    );

    autoPtr<virtual_meshWallInfo>& getVMPlaneInfo
    (
        label ID
    );

    wallContactVars& getWallCntVars()
    {
        return wallCntvar_;
    }

    void clearWallCntVars();

    scalar getMu(wallContactVars& wallCntvar)
    {
        return wallCntvar.getMeanCntPar().aMu_;
    }

    vector getLVec
    (
        wallContactVars& wallCntvar,
        ibContactClass ibCClass
    );

    vector getVeli
    (
        wallContactVars& wallCntvar,
        ibContactVars& cVars
    );

    void evalVariables
    (
        wallContactVars& wallCntvar,
        ibContactClass& ibCClass,
        ibContactVars& cVars
    );

    vector getFNe(wallContactVars& wallCntvar);

    vector getFA(wallContactVars& wallCntvar);

    vector getFNd(wallContactVars& wallCntvar);

    vector getFt(wallContactVars& wallCntvar, scalar deltaT);

    forces& get_out_force()
    {
        return outForce_;
    }

    void setResolvedContact(bool resolved)
    {
        contact_resolved_ = resolved;
    }

    bool getContactResolved()
    {
        return contact_resolved_;
    }

    bool isContactPatch(string patchName)
    {
        return wallMeanPars_.found(patchName);
    }

    HashTable<physicalProperties,string,Hash<string>>& getWallMeanPars()
    {
        return wallMeanPars_;
    }

    List<autoPtr<virtual_meshWallInfo>>& getVMContactInfo()
    {
        return vmWInfoList_;
    }

    List<autoPtr<virtual_meshWallInfo>>& getVMPlaneInfo()
    {
        return vmPlaneInfoList_;
    }

    const List<Tuple2<point,boundBox>>& getInternalElements()
    {
        return internalBBData_;
    }

    const List<string>& getContactPatches()
    {
        return contactPatches_;
    }

    void syncData();

    void syncContactResolve();

    label& getBodyId()
    {
        return body_id_;
    }

    boundBox& getsWCBB()
    {
        return BB_;
    }

};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //