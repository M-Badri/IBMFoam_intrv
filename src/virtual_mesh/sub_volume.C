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
#include "sub_volume.H"

using namespace Foam;

//---------------------------------------------------------------------------//
sub_volume::sub_volume()
{
}

sub_volume::sub_volume
(
    const boundBox bb
)
:
treeBoundBox(bb),
parentSV_(nullptr),
cVolumeInfo_(volumeType::unknown),
tVolumeInfo_(volumeType::unknown)
{
}

sub_volume::sub_volume
(
    const boundBox bb,
    const std::shared_ptr<sub_volume> parentSV,
    const volumeType cVolumeType,
    const volumeType tVolumeType
)
:
treeBoundBox(bb),
parentSV_(parentSV),
cVolumeInfo_(cVolumeType),
tVolumeInfo_(tVolumeType),
isEdge_(false)
{
}

sub_volume::~sub_volume()
{
}
//---------------------------------------------------------------------------//
ibSubVolumeInfo& sub_volume::cVolumeInfo()
{
    return cVolumeInfo_;
}
//---------------------------------------------------------------------------//
ibSubVolumeInfo& sub_volume::tVolumeInfo()
{
    return tVolumeInfo_;
}
//---------------------------------------------------------------------------//
List<sub_volume>& sub_volume::childSubVolumes()
{
    if (childSubVolumes_.size() == 0)
    {
        std::shared_ptr<sub_volume> parentSV = std::make_shared<sub_volume>(*this);
        for (direction octant = 0; octant < 8; octant++)
        {
            childSubVolumes_.append
            (
                sub_volume
                (
                    subBbox(octant),
                    parentSV,
                    cVolumeInfo_.volumeType_,
                    tVolumeInfo_.volumeType_
                )
            );
        }
    }

    return childSubVolumes_;
}
//---------------------------------------------------------------------------//
bool sub_volume::hasChildSubVolumes() const
{
    return childSubVolumes_.size() > 0;
}
//---------------------------------------------------------------------------//
std::shared_ptr<sub_volume>& sub_volume::parentSV()
{
    return parentSV_;
}
//---------------------------------------------------------------------------//
ibSubVolumeInfo& sub_volume::getVolumeInfo(bool cIb)
{
    if (cIb)
    {
        return cVolumeInfo_;
    }
    else
    {
        return tVolumeInfo_;
    }
}
//---------------------------------------------------------------------------//
const ibSubVolumeInfo& sub_volume::getVolumeInfo(bool cIb) const
{
    if (cIb)
    {
        return cVolumeInfo_;
    }
    else
    {
        return tVolumeInfo_;
    }
}
//---------------------------------------------------------------------------//
void sub_volume::setAsEdge()
{
    isEdge_ = true;
}
//---------------------------------------------------------------------------//
bool sub_volume::isEdge() const
{
    return isEdge_;
}
//---------------------------------------------------------------------------//
