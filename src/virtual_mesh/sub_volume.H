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
\*---------------------------------------------------------------------------*/

#ifndef sub_volume_H
#define sub_volume_H

#include "treeBoundBox.H"
#include "volumeType.H"

#include <memory>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

struct ibSubVolumeInfo
{
    ibSubVolumeInfo():
        volumeType_(volumeType::unknown)
    {}

    ibSubVolumeInfo(volumeType vType)
    {
        switch(vType)
        {
            case volumeType::unknown:
            case volumeType::mixed:
                volumeType_ = volumeType::unknown;
                break;
            case volumeType::inside:
                volumeType_ = volumeType::inside;
                break;
            case volumeType::outside:
                volumeType_ = volumeType::outside;
                break;
        }
    }

    volumeType volumeType_;
    autoPtr<labelList> shapesIn_;
};

/*---------------------------------------------------------------------------*\
                        Class sub_volume Declaration
\*---------------------------------------------------------------------------*/
class sub_volume : public treeBoundBox
{
    private:

    std::shared_ptr<sub_volume> parentSV_;
    List<sub_volume> childSubVolumes_;

    ibSubVolumeInfo cVolumeInfo_;
    ibSubVolumeInfo tVolumeInfo_;

    bool isEdge_;

    public:

    // Constructors ------------------------------------------------
    sub_volume();

    sub_volume
    (
        const boundBox bb
    );

    sub_volume
    (
        const boundBox bb,
        const std::shared_ptr<sub_volume> parentSV,
        const volumeType cVolumeType,
        const volumeType tVolumeType
    );

    // Destructors -------------------------------------------------
    ~sub_volume();

    // Member Functions --------------------------------------------

    ibSubVolumeInfo& cVolumeInfo();
    ibSubVolumeInfo& tVolumeInfo();

    List<sub_volume>& childSubVolumes();
    bool hasChildSubVolumes() const;

    std::shared_ptr<sub_volume>& parentSV();

    ibSubVolumeInfo& getVolumeInfo(bool cIb);
    const ibSubVolumeInfo& getVolumeInfo(bool cIb) const;

    void setAsEdge();
    bool isEdge() const;
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
