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
    class for storage of contact variables for wall

SourceFiles
    material_info.C
\*---------------------------------------------------------------------------*/

#ifndef material_info_H
#define material_info_H

#include "dictionary.H"
#include "fvCFD.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                          Class material_info Declaration
\*---------------------------------------------------------------------------*/

class material_info
{
    private:

    typedef DynamicList<label> DynamicLabelList;
    typedef DynamicList<scalar> DynamicScalarList;
    typedef DynamicList<vector> DynamicVectorList;
    typedef DynamicList<point> DynamicPointList;

    protected:

        string material_;

        
        scalar  Y_;

        
        scalar  nu_;

        scalar mu_;

        
        scalar adh_coef_;

        
        scalar eps_;

    public:

    // Constructors ------------------------------------------------
    material_info
    (
        string material,
        scalar Y,
        scalar nu,
        scalar mu,
        scalar adh_coef,
        scalar eps
    );

    material_info
    (
        const material_info& mi
    );

    // Destructors -------------------------------------------------
    ~material_info();

    const string& getMaterial() const
    {
        return material_;
    }

    const scalar& getY()   const
    {
        return Y_;
    }

    const scalar& getNu()   const
    {
        return nu_;
    }

    const scalar& getMu()   const
    {
        return mu_;
    }

    const scalar& getAdhN()   const
    {
        return adh_coef_;
    }

    const scalar& getEps()   const
    {
        return eps_;
    }
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
