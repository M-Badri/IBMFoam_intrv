/*---------------------------------------------------------------------------*\
License

    IBMFoam is licensed under the GNU LESSER GENERAL PUBLIC LICENSE (LGPL).

    Everyone is permitted to copy and distribute verbatim copies of this license
    document, but changing it is not allowed.

    This version of the GNU Lesser General Public License incorporates the terms
    and conditions of version 3 of the GNU General Public License, supplemented
    by the additional permissions listed below.

    You should have received a copy of the GNU Lesser General Public License
    along with openHFDIB. If not, see <   <http://www.gnu.org/licenses/lgpl.html>>.

InNamspace
    Foam
\*---------------------------------------------------------------------------*/
#include "virtual_meshWall.H"

#include "virtual_meshLevel.H"

using namespace Foam;

//---------------------------------------------------------------------------//
virtual_meshWall::virtual_meshWall
(
    virtual_meshWallInfo& vMeshWallInfo,
    geom_model& cGeomModel
)
:
cGeomModel_(cGeomModel),
vMeshWallInfo_(vMeshWallInfo),
bb_matrix_(vMeshWallInfo.sub_volumeNVector,
    vMeshWallInfo.bBox,
    vMeshWallInfo.charCellSize,
    vMeshWallInfo.sub_volumeV)
{}

virtual_meshWall::~virtual_meshWall()
{
}
//---------------------------------------------------------------------------//
bool virtual_meshWall::detectFirstContactPoint()
{
    autoPtr<DynamicVectorList> nextToCheck(
        new DynamicVectorList);

    autoPtr<DynamicVectorList> auxToCheck(
        new DynamicVectorList);

    nextToCheck->append(bb_matrix_.getSVIndexForPoint_Wall(vMeshWallInfo_.getStartingPoint()));
    nextToCheck->append(bb_matrix_.cornerNeighbourSubVolumes(nextToCheck()[0]));
    // InfoH << DEM_Info << " -- VM firstSV : " << nextToCheck()[0] << " point " << bb_matrix_[nextToCheck()[0]].center << endl;
    while (nextToCheck->size() > 0)
    {
        auxToCheck->clear();
        forAll (nextToCheck(),sV)
        {
            sub_volumeProperties& cSubVolume = bb_matrix_[nextToCheck()[sV]];
            if (!cSubVolume.toCheck)
            {
                continue;
            }
            checkSubVolume(cSubVolume);

            if (cSubVolume.isCBody)
            {
                vMeshWallInfo_.startingPoint = cSubVolume.center;
                resetSubVolume(cSubVolume);

                return true;

            }
            auxToCheck().append(bb_matrix_.cornerNeighbourSubVolumes(nextToCheck()[sV]));
        }
        const autoPtr<DynamicVectorList> helpPtr(nextToCheck.ptr());
        nextToCheck.set(auxToCheck.ptr());
        auxToCheck = helpPtr;
    }
    return false;
}
//---------------------------------------------------------------------------////---------------------------------------------------------------------------//
bool virtual_meshWall::detectFirstFaceContactPoint()
{
    autoPtr<DynamicVectorList> nextToCheck(
        new DynamicVectorList);

    autoPtr<DynamicVectorList> auxToCheck(
        new DynamicVectorList);

    nextToCheck->append(bb_matrix_.getSVIndexForPoint_Wall(vMeshWallInfo_.getStartingPoint()));
    nextToCheck->append(bb_matrix_.faceNeighbourSubVolumes(nextToCheck()[0]));
    // InfoH << DEM_Info << " -- VM firstSV : " << nextToCheck()[0] << " point " << bb_matrix_[nextToCheck()[0]].center << endl;
    while (nextToCheck->size() > 0)
    {
        auxToCheck->clear();
        forAll (nextToCheck(),sV)
        {
            sub_volumeProperties& cSubVolume = bb_matrix_[nextToCheck()[sV]];
            if (!cSubVolume.toCheck)
            {
                continue;
            }
            checkSubVolume(cSubVolume);

            if (cSubVolume.isCBody)
            {
                vMeshWallInfo_.startingPoint = cSubVolume.center;
                resetSubVolume(cSubVolume);

                return true;

            }
            auxToCheck().append(bb_matrix_.faceNeighbourSubVolumes(nextToCheck()[sV]));
        }
        const autoPtr<DynamicVectorList> helpPtr(nextToCheck.ptr());
        nextToCheck.set(auxToCheck.ptr());
        auxToCheck = helpPtr;
    }
    return false;
}
//---------------------------------------------------------------------------//
scalar virtual_meshWall::evaluateContact()
{
    label volumeCount = 0;
    contactCenter_ = vector::zero;
    autoPtr<DynamicVectorList> nextToCheck(
        new DynamicVectorList);
    autoPtr<DynamicVectorList> auxToCheck(
        new DynamicVectorList);
    nextToCheck->append(bb_matrix_.getSVIndexForPoint_Wall(vMeshWallInfo_.getStartingPoint()));
    label iterCount(0);
    while (nextToCheck().size() > 0)
    {
        auxToCheck().clear();

        forAll (nextToCheck(),sV)
        {
            sub_volumeProperties& cSubVolume = bb_matrix_[nextToCheck()[sV]];
            iterCount++;
            if (!cSubVolume.toCheck)
            {
                continue;
            }

            checkSubVolume(cSubVolume);
            if (cSubVolume.isCBody)
            {
                volumeCount++;
                contactCenter_ += cSubVolume.center;
                auxToCheck->append(bb_matrix_.faceNeighbourSubVolumes(nextToCheck()[sV]));
            }
        }
        const autoPtr<DynamicVectorList> helpPtr(nextToCheck.ptr());
        nextToCheck.set(auxToCheck.ptr());
        auxToCheck = helpPtr;
    }
    if (volumeCount > 0)
    {
        contactCenter_ /= volumeCount;
    }

    return volumeCount*bb_matrix_.getSubVolumeV();
}
//---------------------------------------------------------------------------//
void virtual_meshWall::checkSubVolume(sub_volumeProperties& sub_volume)
{
    if (sub_volume.toCheck)
    {
        sub_volume.isCBody = cGeomModel_.pointInside(sub_volume.center);
        sub_volume.toCheck = false;
    }
}
//---------------------------------------------------------------------------//
void virtual_meshWall::resetSubVolume(sub_volumeProperties& sub_volume)
{
    sub_volume.toCheck = true;
    sub_volume.isCBody = false;
    sub_volume.isTBody = false;
    sub_volume.isOnEdge = false;   
}
//---------------------------------------------------------------------------//
label virtual_meshWall::getInternalSV()
{
    autoPtr<DynamicVectorList> nextToCheck(
        new DynamicVectorList);
    autoPtr<DynamicVectorList> auxToCheck(
        new DynamicVectorList);

    nextToCheck->append(bb_matrix_.getSVIndexForPoint(vMeshWallInfo_.getStartingPoint()));
    label innerSVCount(0);
    vectorHashSet octreeSvSet;
    
    while (nextToCheck->size() > 0)
    {
        auxToCheck().clear();
        forAll (nextToCheck(),sV)
        {   
            if (!octreeSvSet.found(nextToCheck()[sV]))
            {
                octreeSvSet.insert(nextToCheck()[sV]);
                sub_volumeProperties& cSubVolume = bb_matrix_[nextToCheck()[sV]];
                if (cSubVolume.isCBody)
                {
                    bool isNotOnEdge(true);
                    List<vector> neighbourSubVolumes = bb_matrix_.faceNeighbourSubVolumes(nextToCheck()[sV]);
                    neighbourSubVolumes.append(bb_matrix_.edgeNeighbourSubVolumes(nextToCheck()[sV]));
                    neighbourSubVolumes.append(bb_matrix_.cornerNeighbourSubVolumes(nextToCheck()[sV]));
                    forAll(neighbourSubVolumes,nSV)
                    {
                        isNotOnEdge *= bb_matrix_[neighbourSubVolumes[nSV]].isCBody;
                    }
                    if(!isNotOnEdge)
                    {
                        cSubVolume.isOnEdge = true;
                    }
                    auxToCheck().append(bb_matrix_.faceNeighbourSubVolumes(nextToCheck()[sV]));
                }
                if(!cSubVolume.isOnEdge)
                {
                    innerSVCount++;
                }
            }
        }

        const autoPtr<DynamicVectorList> helpPtr(nextToCheck.ptr());
        nextToCheck.set(auxToCheck.ptr());
        auxToCheck = helpPtr;
    }
    return innerSVCount;
}
// ************************************************************************* //
