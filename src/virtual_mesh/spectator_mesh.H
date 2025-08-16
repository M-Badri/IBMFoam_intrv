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
    wallContactInfo.C

\*---------------------------------------------------------------------------*/

#ifndef spectator_mesh_H
#define spectator_mesh_H

#include "dictionary.H"
#include "fvCFD.H"
#include "ibContactClass.H"
#include "sMProperties.H"
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                        Class ibContactClass Declaration
\*---------------------------------------------------------------------------*/
class spectator_mesh
{
    private:

    typedef DynamicList<label> DynamicLabelList;
    typedef DynamicList<scalar> DynamicScalarList;
    typedef DynamicList<vector> DynamicVectorList;
    typedef DynamicList<point> DynamicPointList;

    public:
    
    const vector matrixSize_;
    const boundBox bBox_;
    const scalar charCellSize_;

    List<List<List<autoPtr<sMProperties>>>> centroidMatrix_;
    List<List<List<autoPtr<sMProperties>>>> vertexMatrix_;
    // Constructors ------------------------------------------------

    spectator_mesh
    (
    const vector matrixSize,
    const boundBox bBox,
    const scalar charCellSize
    );

    // Destructors -------------------------------------------------
    ~spectator_mesh();

//-----------------------------operator-------------------------------//
    sMProperties& operator[](const vector& sVI)
    {
        if (!centroidMatrix_[sVI[0]][sVI[1]][sVI[2]].valid())
        {
            centroidMatrix_[sVI[0]][sVI[1]][sVI[2]].set(
                new sMProperties(getCentroidPoint(sVI))
            );
        }
        return centroidMatrix_[sVI[0]][sVI[1]][sVI[2]]();
    }
//-----------------------------operator-------------------------------// 
   sMProperties& operator()(const vector& sVI)
    {
        if (!vertexMatrix_[sVI[0]][sVI[1]][sVI[2]].valid())
        {
            vertexMatrix_[sVI[0]][sVI[1]][sVI[2]].set(
                new sMProperties(getVertexPoint(sVI))
            );
        }
        return vertexMatrix_[sVI[0]][sVI[1]][sVI[2]]();
    }
//-----------------------------operator-------------------------------//
//-----------------------------Func-------------------------------//
    point getCentroidPoint
    (
        const vector& elementIndex
    );

    point getVertexPoint
    (
        const point& elementIndex
    );

    vector getSMCentroidIndex
    (
        const point& pointInDomain
    );

    vector getSMVertexIndex
    (
        const point& pointInDomain
    );

    List<vector> faceNeighbourElements
    (
        vector& elementIndex
    );    

    List<vector> elementVertexIndexies
    (
        vector& elementIndex
    );

    bool isPointInElementBB
    (
        point& pointToCheck,
        vector& elementIndex
    );

    boundBox getElementBB
    (
        vector& elementIndex
    );

};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
