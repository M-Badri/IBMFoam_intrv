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
    class for Verlet point
\*---------------------------------------------------------------------------*/

#ifndef verletPoint_H
#define verletPoint_H

#include "dictionary.H"
#include "IDLList.H"
#include "fvCFD.H"
#include "geom_model.H"

#include <memory>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{


class verletBox;

class verletPoint
{
    private:

    label bodyID_;
    std::shared_ptr<verletBox> parentBox_;

    point& vPoint_;
    const bool isMin_;
    bool isStatic_;

    public:

    // Constructors ------------------------------------------------
    verletPoint
    (
        label bodyID,
        std::shared_ptr<verletBox> parentBox,
        point& vPoint,
        const bool isMin,
        bool isStatic
    )
    :
    bodyID_(bodyID),
    parentBox_(parentBox),
    vPoint_(vPoint),
    isMin_(isMin),
    isStatic_(isStatic)
    {};

    // Destructors -------------------------------------------------
    virtual ~verletPoint(){};

    bool isMin()
    {
        return isMin_;
    }

    point& getPoint()
    {
        return vPoint_;
    }

    label getBodyId()
    {
        return bodyID_;
    }

    bool getIsStatic()
    {
        return isStatic_;
    }

    std::shared_ptr<verletBox> getParentBox()
    {
        return parentBox_;
    }
};

class verletBox
: public std::enable_shared_from_this<verletBox>
{
    private:

    label bodyID_;
    std::shared_ptr<boundBox> bBox_;
    bool isStatic_;

    std::shared_ptr<verletPoint> minPoint_;
    std::shared_ptr<verletPoint> maxPoint_;

    // Constructors ------------------------------------------------
    verletBox
    (
        label bodyID,
        std::shared_ptr<boundBox> bBox,
        bool isStatic
    )
    :
    bodyID_(bodyID),
    bBox_(bBox),
    isStatic_(isStatic)
    {};

    public:

    std::shared_ptr<verletBox> getptr()
    {
        return shared_from_this();
    }
   
    [[nodiscard]] static std::shared_ptr<verletBox> create
    (
        label bodyID,
        std::shared_ptr<boundBox> bBox,
        bool isStatic
    )
    {
        return std::shared_ptr<verletBox>(new verletBox(bodyID, bBox, isStatic));
    }

    std::shared_ptr<boundBox>& getBBox()
    {
        return bBox_;
    }

    label getBodyId()
    {
        return bodyID_;
    }

    void setVerletPoints()
    {
        minPoint_ = std::make_shared<verletPoint>(bodyID_, shared_from_this(), bBox_->min(), true,isStatic_);
        maxPoint_ = std::make_shared<verletPoint>(bodyID_, shared_from_this(), bBox_->max(), false,isStatic_);
    }

    std::shared_ptr<verletPoint> getMinPoint()
    {
        return minPoint_;
    }

    std::shared_ptr<verletPoint> getMaxPoint()
    {
        return maxPoint_;
    }
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
