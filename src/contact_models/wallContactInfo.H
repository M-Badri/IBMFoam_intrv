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
    wallContactInfo.C
\*---------------------------------------------------------------------------*/

#ifndef wallContactInfo_H
#define wallContactInfo_H

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
#include "wallSubContactInfo.H"
#include "wall_plane_info.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                      Class wallContactInfo Declaration
\*---------------------------------------------------------------------------*/

class wallContactInfo
{
    private:

    typedef DynamicList<label> DynamicLabelList;
    typedef DynamicList<scalar> DynamicScalarList;
    typedef DynamicList<vector> DynamicVectorList;
    typedef DynamicList<point> DynamicPointList;
    typedef HashSet<vector,Hash<vector>> vectorHashSet;

    ibContactClass& ibContactClass_;
    ibContactVars& ibContactVars_;

    label body_id_;

    // Mean model parameters
    HashTable<physicalProperties,string,Hash<string>> wallMeanPars_;

    std::vector<std::shared_ptr<wallSubContactInfo>> subCList_;

    scalar reduceM_;   

    forces outForce_; 

    autoPtr<spectator_mesh> SM_;

    autoPtr<virtual_meshWallInfo> vmwInfo_;

    List<string> contactPatches_;

    public:

    // Constructors ------------------------------------------------
    wallContactInfo
    (
        ibContactClass& cClass,
        ibContactVars& cVars
    );

    // Destructors -------------------------------------------------
    ~wallContactInfo();

    HashTable<physicalProperties,string,Hash<string>>& getWallMeanPars()
    {
        return wallMeanPars_;
    }

    ibContactClass& getcClass()
    {
        return ibContactClass_;
    }

    ibContactVars& getcVars()
    {
        return ibContactVars_;
    }

    label getBodyId()
    {
        return body_id_;
    }

    boundBox constructBoundBox
    (
        boundBox& bodyBB
    );

    void constructSM();

    bool isInsidePlane
    (
        vector checkedPoint,
        const string& wall,
        const HashTable<List<vector>,string,Hash<string>>& wallPatches = wall_plane_info::getWallPlaneInfo()
    );

    bool detectWallContact(
        const HashTable<List<vector>,string,Hash<string>>& wallPatches = wall_plane_info::getWallPlaneInfo()
    );

    void findContactAreas();

    List<DynamicList<vector>> detectPossibleSMContact
    (
        DynamicVectorList& contactPoints,
        List<string>& contactPatches
    );

    boundBox getSCBBox
    (
        DynamicVectorList& sub_contactAreas
    );

    boundBox constructVMBox
    (
        boundBox& baseContactAreaBB,
        string& wallName
    );

    vector getPlanePoint
    (
        vector pointInDomain,
        string wallName
    );

    void checkSMElement
    (
        vector& index,
        List<string> contactPatches
    );

    void checkElement
    (
        vector& index,
        bool& inMesh,
        bool& inBody,
        bool& allInMesh,
        bool& allInBody
    );

    boundBox contactPlaneBBox
    (
        boundBox contactBoundBox,
        string wallName
    );

    boundBox correctSMBBforWall
    (
        boundBox bB,
        string wallName
    );

    void clearOldContact();

    void setNewSubContact(
        const List<Tuple2<point,boundBox>>& contactBBData,
        const List<Tuple2<point,boundBox>>& planeBBData,
        const List<string>& contactPatches,
        const List<Tuple2<point,boundBox>>& isInternal,
        boundBox BB
    );

    std::vector<std::shared_ptr<wallSubContactInfo>>& getWallSCList()
    {
        return subCList_;
    }

    void registerSubContactList(DynamicList<wallSubContactInfo*>& contactList);

    List<string>& getContactPatches()
    {
        return contactPatches_;
    }


};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
