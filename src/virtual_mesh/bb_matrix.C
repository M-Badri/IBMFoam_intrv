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
#include "bb_matrix.H"

#include "virtual_meshLevel.H"

using namespace Foam;

//---------------------------------------------------------------------------//
bb_matrix::bb_matrix
(
    const vector matrixSize,
    const boundBox bBox,
    const scalar& charCellSize,
    const scalar& sub_volumeV
)
:
matrixSize_(matrixSize),
bBox_(bBox),
charCellSize_(charCellSize),
sub_volumeV_(sub_volumeV)
{
    bb_matrix_ = List<List<List<autoPtr<sub_volumeProperties>>>>(matrixSize_[0],
        List<List<autoPtr<sub_volumeProperties>>>(matrixSize_[1],
        List<autoPtr<sub_volumeProperties>>(matrixSize_[2])));
}
bb_matrix::~bb_matrix()
{
}
//---------------------------------------------------------------------------//
vector bb_matrix::getPointInMesh
(
	const vector& sub_volumeIndex
)
{
    return vector(
        bBox_.min()[0] + (2*sub_volumeIndex[0]+1)*(charCellSize_/virtual_meshLevel::getLevelOfDivision())*0.5,
        bBox_.min()[1] + (2*sub_volumeIndex[1]+1)*(charCellSize_/virtual_meshLevel::getLevelOfDivision())*0.5,
        bBox_.min()[2] + (2*sub_volumeIndex[2]+1)*(charCellSize_/virtual_meshLevel::getLevelOfDivision())*0.5
    );
}
//---------------------------------------------------------------------------//
vector bb_matrix::getSVIndexForPoint
(
    const point& pointInDomain
)
{
    vector sub_volumeIndex(vector::zero);
    for(label i = 0;i<3;i++)
    {
        sub_volumeIndex[i] = floor((pointInDomain[i]-bBox_.min()[i])/charCellSize_*virtual_meshLevel::getLevelOfDivision());
        if(sub_volumeIndex[i] >= matrixSize_[i] || sub_volumeIndex[i] < 0 )
        {
            sub_volumeIndex = (vector::one)*(-1);
            break;
        }
    }
    return sub_volumeIndex;
}
//---------------------------------------------------------------------------//
vector bb_matrix::getSVIndexForPoint_Wall
(
    point pointInDomain
)
{
    for(label i = 0;i<3;i++)
    {
        if(mag(pointInDomain[i])<SMALL)
        {
            pointInDomain[i] = 0.0;
        }
    } 
    vector sub_volumeIndex(vector::zero);
    for(label i = 0;i<3;i++)
    {
        sub_volumeIndex[i] = floor((pointInDomain[i]-bBox_.min()[i])/charCellSize_*virtual_meshLevel::getLevelOfDivision());
        if(sub_volumeIndex[i] >= matrixSize_[i] || sub_volumeIndex[i] < 0 )
        {
            return(getSVIndexForPoint(bBox_.midpoint()));
            break;
        }
    }
    return sub_volumeIndex;
}
//---------------------------------------------------------------------------//
vector bb_matrix::getFirstSubVolumeIndex
(
    point& sub_volumePoint,
    bool& isInMatrix
)
{
    vector sub_volumeIndex(vector::zero);
    isInMatrix = true;
    for(label i = 0;i<3;i++)
    {
        sub_volumeIndex[i] = floor((sub_volumePoint[i]-bBox_.min()[i])
            /charCellSize_*virtual_meshLevel::getLevelOfDivision());

        if(sub_volumeIndex[i] >= matrixSize_[i] || sub_volumeIndex[i] < 0 )
        {
            sub_volumeIndex = (vector::one)*(-1);
            isInMatrix = false;
            break;
        }
    }
    return sub_volumeIndex;
}
//---------------------------------------------------------------------------//
List<vector> bb_matrix::faceNeighbourSubVolumes
(
    vector& sub_volumeIndex
)
{
    const label i(sub_volumeIndex[0]);
    const label j(sub_volumeIndex[1]);
    const label k(sub_volumeIndex[2]);

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
List<vector> bb_matrix::edgeNeighbourSubVolumes //returns only edgeNeighbours without faceNeighbours
(
    vector& sub_volumeIndex
)
{
    const label i(sub_volumeIndex[0]);
    const label j(sub_volumeIndex[1]);
    const label k(sub_volumeIndex[2]);

    vector edgeNeighbour1(i-1,j,k-1);
    vector edgeNeighbour2(i+1,j,k-1);
    vector edgeNeighbour3(i+1,j,k+1);
    vector edgeNeighbour4(i-1,j,k+1);

    vector edgeNeighbour5(i-1,j+1,k);
    vector edgeNeighbour6(i-1,j-1,k);
    vector edgeNeighbour7(i+1,j-1,k);
    vector edgeNeighbour8(i+1,j+1,k);

    vector edgeNeighbour9(i,j+1,k+1);
    vector edgeNeighbour10(i,j+1,k-1);
    vector edgeNeighbour11(i,j-1,k-1);
    vector edgeNeighbour12(i,j-1,k+1);

    const List<vector> NLL{{edgeNeighbour1,edgeNeighbour2,edgeNeighbour3,edgeNeighbour4,
                            edgeNeighbour5,edgeNeighbour6,edgeNeighbour7,edgeNeighbour8,
                            edgeNeighbour9,edgeNeighbour10,edgeNeighbour11,edgeNeighbour12}};
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
List<vector> bb_matrix::cornerNeighbourSubVolumes
(
    vector& sub_volumeIndex
)
{
    label i(sub_volumeIndex[0]);
    label j(sub_volumeIndex[1]);
    label k(sub_volumeIndex[2]);

    vector cornerNeighbour1(i+1,j+1,k+1);
    vector cornerNeighbour2(i+1,j-1,k+1);
    vector cornerNeighbour3(i-1,j-1,k+1);
    vector cornerNeighbour4(i-1,j+1,k+1);

    vector cornerNeighbour5(i+1,j+1,k-1);
    vector cornerNeighbour6(i+1,j-1,k-1);
    vector cornerNeighbour7(i-1,j-1,k-1);
    vector cornerNeighbour8(i-1,j+1,k-1);

    List<vector> NLL{{cornerNeighbour1,cornerNeighbour2,cornerNeighbour3,cornerNeighbour4,
                    cornerNeighbour5,cornerNeighbour6,cornerNeighbour7,cornerNeighbour8}};

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
