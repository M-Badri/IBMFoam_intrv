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
    prtContact.C
\*---------------------------------------------------------------------------*/

#ifndef prtPrtContact_H
#define prtPrtContact_H

#include "dictionary.H"
#include "fvCFD.H"


#include "ibContactClass.H"
#include "prtContactVars.H"
#include "prtContactInfo.H"
#include "periodicBody.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                 Namespace contactModel functions Declaration
\*---------------------------------------------------------------------------*/

namespace contactModel
{
    typedef DynamicList<label>  DynamicLabelList;
    typedef DynamicList<scalar> DynamicScalarList;
    typedef DynamicList<vector> DynamicVectorList;
    typedef DynamicList<point>  DynamicPointList;

    void getContacts(
        const fvMesh&   mesh,
        prtContactInfo& prtcInfo
    );

    void getContacts_Cluster(
        const fvMesh&   mesh,
        prtContactInfo& prtcInfo
    );

    bool detectPrtPrtContact
    (
        const fvMesh& mesh,
        ibContactClass& cClass,
        ibContactClass& tClass,
        prtSubContactInfo& subCInfo
    );

    void getContacts_ArbShape(
        const fvMesh&   mesh,
        prtContactInfo& prtcInfo
    );

    bool detectPrtPrtContact_ArbShape
    (
        const fvMesh& mesh,
        ibContactClass& cClass,
        ibContactClass& tClass,
        prtSubContactInfo& subCInfo
    );

    bool detectPrtPrtContact_Sphere
    (
        const fvMesh& mesh,
        ibContactClass& cClass,
        ibContactClass& tClass
    );

    bool detectPrtPrtContact_Cluster
    (
        const fvMesh& mesh,
        ibContactClass& cClass,
        ibContactClass& tClass,
        prtSubContactInfo& subCInfo
    );

    void getPrtContactVars(
        const fvMesh&   mesh,
        ibContactClass& cClass,
        ibContactClass& tClass,
        prtSubContactInfo& subCInfo
    );

    void getPrtContactVars_ArbShape(
        const fvMesh&   mesh,
        ibContactClass& cClass,
        ibContactClass& tClass,
        prtSubContactInfo& subCInfo
    );

    void getPrtContactVars_Sphere
    (
        const fvMesh&   mesh,
        ibContactClass& cClass,
        ibContactClass& tClass,
        prtSubContactInfo& subCInfo
    );

    void getPrtContactVars_Cluster
    (
        const fvMesh&   mesh,
        ibContactClass& cClass,
        ibContactClass& tClass,
        prtSubContactInfo& subCInfo
    );

    bool solvePrtContact(
        const fvMesh&   mesh,
        prtContactInfo& cInfo,
        prtSubContactInfo& subCInfo,
        scalar deltaT
    );
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
