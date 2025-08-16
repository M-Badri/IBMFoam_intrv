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

\*---------------------------------------------------------------------------*/
#include "spectator_mesh.H"

using namespace Foam;

//---------------------------------------------------------------------------//
spectator_mesh::spectator_mesh
(
    const vector matrixSize,
    const boundBox bBox,
    const scalar charCellSize
)
:
matrixSize_(matrixSize),
bBox_(bBox),
charCellSize_(charCellSize)
{
    // InfoH << DEM_Info << "SM is Alive" << endl;
    centroidMatrix_ = List<List<List<autoPtr<sMProperties>>>>(matrixSize[0],
        List<List<autoPtr<sMProperties>>>(matrixSize[1],
        List<autoPtr<sMProperties>>(matrixSize[2])));

    vertexMatrix_ = List<List<List<autoPtr<sMProperties>>>>(matrixSize[0]+1,
        List<List<autoPtr<sMProperties>>>(matrixSize[1]+1,
        List<autoPtr<sMProperties>>(matrixSize[2]+1)));
    // InfoH << DEM_Info << "SM is Alive I " << endl;
}
spectator_mesh::~spectator_mesh()
{
}
//---------------------------------------------------------------------------//
point spectator_mesh::getCentroidPoint
(
	const vector& elementIndex
)
{
    return vector(
        bBox_.min()[0] + (2*elementIndex[0]+1)*(charCellSize_)*0.5,
        bBox_.min()[1] + (2*elementIndex[1]+1)*(charCellSize_)*0.5,
        bBox_.min()[2] + (2*elementIndex[2]+1)*(charCellSize_)*0.5
    );
}
//---------------------------------------------------------------------------//
point spectator_mesh::getVertexPoint
(
	const vector& elementIndex
)
{
    return vector(
        bBox_.min()[0] + (elementIndex[0])*(charCellSize_),
        bBox_.min()[1] + (elementIndex[1])*(charCellSize_),
        bBox_.min()[2] + (elementIndex[2])*(charCellSize_)
    );
}
//---------------------------------------------------------------------------//
vector spectator_mesh::getSMCentroidIndex
(
    const point& pointInDomain
)
{
    vector centroidIndex(vector::zero);
    for(label i = 0;i<3;i++)
    {
        centroidIndex[i] = floor((pointInDomain[i]-bBox_.min()[i])/charCellSize_);
        if(centroidIndex[i] >= matrixSize_[i] || centroidIndex[i] < 0 )
        {
            centroidIndex = (vector::one)*(-1);
            break;
        }
    }
    return centroidIndex;
}
//---------------------------------------------------------------------------//
vector spectator_mesh::getSMVertexIndex
(
    const point& pointInDomain
)
{
    vector vertexIndex(vector::zero);
    for(label i = 0;i<3;i++)
    {
        vertexIndex[i] = floor((pointInDomain[i]-bBox_.min()[i])/charCellSize_);
        if(vertexIndex[i] >= matrixSize_[i] +1 || vertexIndex[i] < 0 )
        {
            vertexIndex = (vector::one)*(-1);
            break;
        }
    }
    return vertexIndex;
}
//---------------------------------------------------------------------------//
List<vector> spectator_mesh::faceNeighbourElements
(
    vector& elementIndex
)
{
    const label i(elementIndex[0]);
    const label j(elementIndex[1]);
    const label k(elementIndex[2]);

    vector neighbour1(i-1,j,k);
    vector neighbour2(i+1,j,k);
    vector neighbour3(i,j-1,k);
    vector neighbour4(i,j+1,k);
    vector neighbour5(i,j,k-1);
    vector neighbour6(i,j,k+1);

    const List<vector> NLL{{neighbour1,neighbour2,neighbour3,neighbour4,neighbour5,neighbour6}};
    List<vector> cellNeighbours;

    forAll(NLL, neig)
    {
        if((NLL[neig][0] >= 0 && NLL[neig][0] < matrixSize_[0])&&
            (NLL[neig][1] >= 0 && NLL[neig][1] < matrixSize_[1])&&
            (NLL[neig][2] >= 0 && NLL[neig][2] < matrixSize_[2]))
        {
            cellNeighbours.append(NLL[neig]);
        }
    }

	return cellNeighbours;
}
//---------------------------------------------------------------------------//
List<vector> spectator_mesh::elementVertexIndexies
(
    vector& elementIndex
)
{
    const label i(elementIndex[0]);
    const label j(elementIndex[1]);
    const label k(elementIndex[2]);

    List<vector> indexies;

    for(int a = 0; a<2; a++)
    {
        for(int b = 0; b<2; b++)
        {
            for(int c = 0; c<2; c++)
            {
                if((a+i >= 0 && a+i < matrixSize_[0]+1)&&
                    (b+j >= 0 && b+j < matrixSize_[1]+1)&&
                    (c+k >= 0 && c+k < matrixSize_[2]+1))
                {
                    indexies.append(vector(a+i,b+j,c+k));
                }
            }
        }
    } 

    return indexies;
}
//---------------------------------------------------------------------------//
bool spectator_mesh::isPointInElementBB
(
    point& pointToCheck,
    vector& elementIndex
)
{
    List<vector> vertexIndex = elementVertexIndexies(elementIndex);
    pointField elementPoints;
    forAll(vertexIndex, vI)
    {
        const vector lVI = vertexIndex[vI];
        elementPoints.append(vertexMatrix_[lVI[0]][lVI[1]][lVI[2]]().center);
    }
    boundBox elementBB = boundBox(elementPoints,false);

    return elementBB.contains(pointToCheck);
}
//---------------------------------------------------------------------------//
boundBox spectator_mesh::getElementBB
(
    vector& elementIndex
)
{
    List<vector> vertexIndex = elementVertexIndexies(elementIndex);
    pointField elementPoints;
    forAll(vertexIndex, vI)
    {
        const vector lVI = vertexIndex[vI];
        elementPoints.append(vertexMatrix_[lVI[0]][lVI[1]][lVI[2]]().center);
    }
    return boundBox(elementPoints,false);
}
//---------------------------------------------------------------------------//