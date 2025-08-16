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

#ifndef virtual_meshWall_H
#define virtual_meshWall_H

#include "dictionary.H"
#include "fvCFD.H"
#include "geom_model.H"
#include "bb_matrix.H"
#include "clockTime.H"
#include "sub_volumeProperties.H"
// #include ""


// #include "fvMesh.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

struct virtual_meshWallInfo
{
    virtual_meshWallInfo():
    bBox(boundBox()),
    startingPoint(vector::zero),
    sub_volumeNVector(vector::zero),
    charCellSize(0),
    sub_volumeV(0)

    {}
    virtual_meshWallInfo(
        boundBox bBox,
        point startingPoint,
        vector sub_volumeNVector,
        scalar charCellSize,
        scalar sub_volumeV):
    bBox(bBox),
    startingPoint(startingPoint),
    sub_volumeNVector(sub_volumeNVector),
    charCellSize(charCellSize),
    sub_volumeV(sub_volumeV)
    {}

    const point& getStartingPoint()
    {
        return startingPoint;
    }

    const scalar& getSVVolume()
    {
        return sub_volumeV;
    }

    const scalar& getCharCellSize()
    {
        return charCellSize;
    }

    void shiftBBox
    (
        vector shiftVec
    )
    {
        bBox.min() += shiftVec;
        bBox.max() += shiftVec;
    }

    boundBox& getbBox()
    {
        return bBox;
    }

    void setStartingPoint
    (
        point sPoint
    )
    {
        startingPoint = sPoint;
    }

    boundBox bBox;
    point startingPoint;
    vector sub_volumeNVector;
    scalar charCellSize;
    scalar sub_volumeV;    
};

/*---------------------------------------------------------------------------*\
                        Class virtual_meshWall Declaration
\*---------------------------------------------------------------------------*/
class virtual_meshWall
{
    private:

    typedef DynamicList<label> DynamicLabelList;
    typedef DynamicList<scalar> DynamicScalarList;
    typedef DynamicList<vector> DynamicVectorList;
    typedef DynamicList<point> DynamicPointList;
    typedef HashSet<vector,Hash<vector>> vectorHashSet;

    geom_model& cGeomModel_;

    virtual_meshWallInfo& vMeshWallInfo_;

    bb_matrix bb_matrix_;

    DynamicPointList edgeSubVolumesPoints_;

    public:
    
    vector contactCenter_ = vector::zero;

    // Constructors ------------------------------------------------

    virtual_meshWall
    (
        virtual_meshWallInfo& vMeshWallInfo,
        geom_model& cGeomModel
    );


    // Destructors -------------------------------------------------
    ~virtual_meshWall();

    bool detectFirstContactPoint();

    bool detectFirstFaceContactPoint();

    scalar evaluateContact();

    void checkSubVolume(sub_volumeProperties& sub_volume);

    void resetSubVolume(sub_volumeProperties& sub_volume);

    vector shiftEdgeSubVolume
    (
        vector& sub_volumeIndex
    );

    label getInternalSV();

    DynamicPointList& getEdgeSVPoints()
    {
        return edgeSubVolumesPoints_;
    }

    vector& getContactCenter()
    {
        return contactCenter_;
    }

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
