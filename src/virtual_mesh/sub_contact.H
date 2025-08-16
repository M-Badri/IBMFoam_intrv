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
    cass storing sub_contact information
--------------------------------------------*/

#ifndef sub_contact_H
#define sub_contact_H

#include "sub_volume.H"

#include <memory>
#include <vector>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                        Class sub_contact Declaration
\*---------------------------------------------------------------------------*/
class sub_contact
{
    private:

    std::vector<std::shared_ptr<sub_volume>> sub_volumes_;

    boundBox boundBox_;
    scalar volume_;

    public:

    // Constructors ------------------------------------------------
    sub_contact();

    // Destructors -------------------------------------------------
    ~sub_contact();

    // Member Functions --------------------------------------------

    void addSubVolume(std::shared_ptr<sub_volume> sV);

    bool canCombine(sub_volume& sV);

    const boundBox& getBounds() const
    {
        return boundBox_;
    }

    scalar getVolume() const
    {
        return volume_;
    }

    const std::vector<std::shared_ptr<sub_volume>>& getSubVolumes() const
    {
        return sub_volumes_;
    }

    DynamicList<point> getEdgePoints() const;
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
