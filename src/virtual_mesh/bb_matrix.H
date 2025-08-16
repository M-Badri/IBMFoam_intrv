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

#ifndef bb_matrix_H
#define bb_matrix_H

#include "dictionary.H"
#include "fvCFD.H"
#include "ibContactClass.H"
#include "sub_volumeProperties.H"
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                        Class ibContactClass Declaration
\*---------------------------------------------------------------------------*/
class bb_matrix
{
    private:

    typedef DynamicList<label> DynamicLabelList;
    typedef DynamicList<scalar> DynamicScalarList;
    typedef DynamicList<vector> DynamicVectorList;
    typedef DynamicList<point> DynamicPointList;

    const vector matrixSize_;
    const boundBox bBox_;
    const scalar charCellSize_;
    const scalar sub_volumeV_;
    List<List<List<autoPtr<sub_volumeProperties>>>> bb_matrix_;

    public:

    // Constructors ------------------------------------------------

    bb_matrix
    (
    const vector matrixSize,
    const boundBox bBox,
    const scalar& charCellSize,
    const scalar& sub_volumeV
    );

    // Destructors -------------------------------------------------
    ~bb_matrix();

//-----------------------------operator-------------------------------//
    sub_volumeProperties& operator[](const vector& sVI)
    {
        if (!bb_matrix_[sVI[0]][sVI[1]][sVI[2]].valid())
        {
            bb_matrix_[sVI[0]][sVI[1]][sVI[2]].set(
                new sub_volumeProperties(getPointInMesh(sVI))
            );
        }

        return bb_matrix_[sVI[0]][sVI[1]][sVI[2]]();
    }
//-----------------------------operator-------------------------------//
//-----------------------------Func-------------------------------//
    vector getPointInMesh
    (
        const vector& sub_volumeIndex
    );

    vector getSVIndexForPoint
    (
        const point& pointInDomain
    );

    vector getSVIndexForPoint_Wall
    (
        point pointInDomain
    );

    vector getFirstSubVolumeIndex
    (
        point& sub_volumePoint,
        bool& isInMatrix
    );

    List<vector> faceNeighbourSubVolumes
    (
        vector& sub_volumeIndex
    );

    List<vector> edgeNeighbourSubVolumes
    (
        vector& sub_volumeIndex
    );

    List<vector> cornerNeighbourSubVolumes
    (
        vector& sub_volumeIndex
    );

    scalar getSubVolumeV()
    {
        return sub_volumeV_;
    }

    scalar getCharCellSize()
    {
        return charCellSize_;
    }

    boundBox getBBox()
    {
        return bBox_;
    }
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
