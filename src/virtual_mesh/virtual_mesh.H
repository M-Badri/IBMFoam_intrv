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

#ifndef virtual_mesh_H
#define virtual_mesh_H

#include "dictionary.H"
#include "fvCFD.H"
#include "geom_model.H"
#include "clockTime.H"
#include "sub_volume.H"
#include "sub_contact.H"


// #include "fvMesh.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

struct virtual_meshInfo
{
    virtual_meshInfo():
    sV(boundBox()),
    sub_volumeV(0)
    {}

    virtual_meshInfo(const virtual_meshInfo& vMeshInfo):
    sV(vMeshInfo.sV),
    sub_volumeV(vMeshInfo.sub_volumeV)
    {
        if (vMeshInfo.startingPoint.valid())
        {
            startingPoint.set(new point(vMeshInfo.startingPoint()));
        }
    }

    virtual_meshInfo(
        boundBox& bBox,
        scalar sub_volumeV):
    sV(bBox),
    sub_volumeV(sub_volumeV)
    {}

    void operator= (const virtual_meshInfo& vMeshInfo)
    {
        sV = vMeshInfo.sV;
        sub_volumeV = vMeshInfo.sub_volumeV;

        if (vMeshInfo.startingPoint.valid())
        {
            startingPoint.reset(new point(vMeshInfo.startingPoint()));
        }
    }

    const point& getStartingPoint()
    {
        if (!startingPoint.valid())
        {
            startingPoint.set(new point(sV.midpoint()));
        }

        return startingPoint();
    }

    sub_volume& getSubVolume()
    {
        return sV;
    }

    sub_volume sV;
    scalar sub_volumeV;
    autoPtr<point> startingPoint;
};

/*---------------------------------------------------------------------------*\
                        Class virtual_mesh Declaration
\*---------------------------------------------------------------------------*/
class virtual_mesh
{
    private:

    typedef DynamicList<label> DynamicLabelList;
    typedef DynamicList<scalar> DynamicScalarList;
    typedef DynamicList<vector> DynamicVectorList;
    typedef DynamicList<point> DynamicPointList;
    typedef HashSet<vector,Hash<vector>> vectorHashSet;

    geom_model& cGeomModel_;
    geom_model& tGeomModel_;

    virtual_meshInfo& vMeshInfo_;

    point contactCenter_;
    DynamicPointList edgeSubVolumesPoints_;

    public:


    // Constructors ------------------------------------------------

    virtual_mesh
    (
        virtual_meshInfo& vMeshInfo,
        geom_model& cGeomModel_,
        geom_model& tGeomModel_
    );


    // Destructors -------------------------------------------------
    ~virtual_mesh();

    bool detectFirstContactPoint();

    // Get contact sub_volume
    bool detectFirstVolumeInContact
    (
        sub_volume& sV,
        bool& startPointFound
    );

    scalar evaluateContact();

    void inspectSubVolume(
        sub_volume& sub_volume,
        scalar& contactVolume,
        vector& contactCenter,
        DynamicPointList& edgePoints
    );

    DynamicPointList& getEdgeSVPoints()
    {
        return edgeSubVolumesPoints_;
    }

    point& getContactCenter()
    {
        vMeshInfo_.startingPoint() = contactCenter_;
        return contactCenter_;
    }

    Tuple2<scalar,vector> get3DcontactNormalAndSurface(bool nonConvex);

    Tuple2<scalar,vector> get3DcontactNormalAndSurface(DynamicPointList edgeSubVolumesPoints);

    std::vector<sub_contact> findsub_contacts(sub_volume& sV);

    bool canCombineSubContacts(
        sub_contact& main,
        sub_contact& comp
    ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
