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
    class for storage of contact variables

SourceFiles
    ibContactClass.C
\*---------------------------------------------------------------------------*/

#ifndef ibContactClass_H
#define ibContactClass_H

#include "dictionary.H"
#include "fvCFD.H"

#include "geom_model.H"
#include "ibContactVars.H"
#include "wallContactVars.H"
#include "material_info.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                        Class ibContactClass Declaration
\*---------------------------------------------------------------------------*/

class ibContactClass
{
    private:

    typedef DynamicList<label> DynamicLabelList;
    typedef DynamicList<scalar> DynamicScalarList;
    typedef DynamicList<vector> DynamicVectorList;
    typedef DynamicList<point> DynamicPointList;

    protected:

    std::shared_ptr<geom_model>& geom_model_;

   
    bool    isInWallContact_;

    bool inContact_with_static_;

    label timeStepsInContWStatic_;

    const material_info& matInfo_;

    public:

    // Constructors ------------------------------------------------
    ibContactClass
    (
        std::shared_ptr<geom_model>& geom_model,
        const string& material
    );

    ibContactClass
    (
        const ibContactClass& other
    );

    // Destructors -------------------------------------------------
    ~ibContactClass();

    geom_model& get_geom_model()
    {
        return *geom_model_;
    }

    std::shared_ptr<geom_model> get_geom_modelPtr()
    {
        return geom_model_;
    }

    const bool& checkWallContact() const
    {
        return isInWallContact_;
    }

    void setWallContact(bool inContact)
    {
        isInWallContact_ = inContact;
    }

    List<DynamicLabelList>& getSurfCells()
    {
        return geom_model_->getSurfaceCellList();
    }

    List<DynamicLabelList>& getIntCells()
    {
        return geom_model_->getInternalCellList();
    }

    const material_info& getMatInfo()
    {
        return matInfo_;
    }

    const bool& checkInContactWithStatic() const
    {
        return inContact_with_static_;
    }

    void inContact_with_static(bool inContact)
    {
        inContact_with_static_ = inContact;
    }

    const label& getTimeStepsInContWStatic() const
    {
        return timeStepsInContWStatic_;
    }

    void setTimeStepsInContWStatic(label setter)
    {
        timeStepsInContWStatic_ = setter;
    }
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
