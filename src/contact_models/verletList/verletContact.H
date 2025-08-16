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

Description
    class for Verlet contact
\*---------------------------------------------------------------------------*/

#ifndef verletContact_H
#define verletContact_H

#include "verletPoint.H"

#include <list>
#include <memory>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class verletContact
{
    private:

    label cBodyID_;
    label tBodyID_;

    std::list<std::pair<std::shared_ptr<verletBox>, std::shared_ptr<verletBox>>> contactBoxes_;

    public:

    // Constructors ------------------------------------------------
    verletContact():
    cBodyID_(-1),
    tBodyID_(-1)
    {};

    verletContact
    (
        std::pair<label, label> cPair
    )
    :
    cBodyID_(cPair.first),
    tBodyID_(cPair.second)
    {};

    // Destructors -------------------------------------------------
    ~verletContact(){};

    void addContactBox
    (
        std::shared_ptr<verletBox> cBox,
        std::shared_ptr<verletBox> tBox
    )
    {
        for (std::pair<std::shared_ptr<verletBox>, std::shared_ptr<verletBox>>& cB : contactBoxes_)
        {
            if (cB.first == cBox && cB.second == tBox)
            {
                return;
            }
        }

        contactBoxes_.push_back(std::make_pair(cBox, tBox));
    }

    void removeContactBox
    (
        std::shared_ptr<verletBox> cBox,
        std::shared_ptr<verletBox> tBox
    )
    {
        contactBoxes_.remove_if(
            [&cBox, &tBox] (const std::pair<std::shared_ptr<verletBox>, std::shared_ptr<verletBox>>& cB)
            { return cB.first == cBox && cB.second == tBox; });
    }

    std::list<std::pair<std::shared_ptr<verletBox>, std::shared_ptr<verletBox>>>& getContactBoxes()
    {
        return contactBoxes_;
    }

}; // End class verletContact

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
