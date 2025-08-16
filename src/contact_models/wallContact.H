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
    wallContact.C
\*---------------------------------------------------------------------------*/

#ifndef wallContact_H
#define wallContact_H

#include "dictionary.H"
#include "fvCFD.H"

#include "ibContactClass.H"
#include "material_info.H"
#include "wallContactInfo.H"
#include "periodicBody.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

namespace contactModel
{
    typedef DynamicList<label>  DynamicLabelList;
    typedef DynamicList<scalar> DynamicScalarList;
    typedef DynamicList<vector> DynamicVectorList;
    typedef DynamicList<point>  DynamicPointList;

    bool detectWallContact
    (
        const fvMesh&   mesh,
        ibContactClass& ibClass,
        wallContactInfo& wallCntInfo
    );

    bool detectWallContact_ArbShape
    (
        const fvMesh&   mesh,
        ibContactClass& ibClass,
        wallContactInfo& wallCntInfo
    );

    bool detectWallContact_Sphere
    (
        const fvMesh&   mesh,
        ibContactClass& ibClass,
        wallContactInfo& wallCntInfo
    );

    bool detectWallContact_Cluster
    (
        const fvMesh&   mesh,
        ibContactClass& ibClass,
        wallContactInfo& wallCntInfo
    );

    void getWallContactVars
    (
        const fvMesh&   mesh,
        wallContactInfo& wallCntInfo,
        const scalar deltaT,
        wallSubContactInfo& sWC
    );

    void getWallContactVars_ArbShape
    (
        const fvMesh&   mesh,
        wallContactInfo& wallCntInfo,
        const scalar deltaT,
        wallSubContactInfo& sWC
    );

    void getWallContactVars_Sphere
    (
        const fvMesh&   mesh,
        wallContactInfo& wallCntInfo,
        const scalar deltaT,
        wallSubContactInfo& sCW
    );

    void getWallContactVars_Cluster
    (
        const fvMesh&   mesh,
        wallContactInfo& wallCntInfo,
        const scalar deltaT,
        wallSubContactInfo& sWC
    );

    DynamicList<Tuple2<label,string>> getContactFacesSphere
    (
        const fvMesh&   mesh,
        wallContactInfo& wallCntInfo,
        label faceLabel,
        DynamicList<Tuple2<label,string>>& checkedFaces
    );


    bool solveWallContact
    (
        const fvMesh&   mesh,
        wallContactInfo& wallCntInfo,
        scalar deltaT,
        wallSubContactInfo& sWC
    );


    scalar getInterVolume_Sphere
    (
        const fvMesh&   mesh,
        ibContactClass& cClass,
        ibContactVars& cVars,
        plane& cPlane
    );

    scalar sphereContactArea
    (
        const fvMesh&   mesh,
        ibContactClass& cClass,
        ibContactVars& cVars,
        plane& cPlane
    );
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
