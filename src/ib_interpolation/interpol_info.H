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
    class for interpolation

SourceFiles
\*---------------------------------------------------------------------------*/

#ifndef interpol_info_H
#define interpol_info_H

#include "dictionary.H"
#include "fvCFD.H"
#include "triSurfaceMesh.H"
#include "geom_model.H"

namespace Foam
{
    class interpol_info
    {
        typedef DynamicList<label> DynamicLabelList;
        typedef DynamicList<scalar> DynamicScalarList;
        typedef DynamicList<vector> DynamicVectorList;
        typedef DynamicList<point> DynamicPointList;

        protected:

        const Foam::fvMesh& mesh_;
        std::shared_ptr<geom_model>& geom_model_;

        autoPtr<List<point>> ibPointsPtr_; 
        autoPtr<List<vector>> ibNormalsPtr_; 

        public:

        // Constructors ------------------------------------------------
        interpol_info
        (
            const Foam::fvMesh& mesh,
            std::shared_ptr<geom_model>& gModel
        );

        // Destructors -------------------------------------------------
        virtual ~interpol_info();

        List<point>& getIbPoints()
        {
            if(ibPointsPtr_.empty())
            {
                setIntpInfo();
            }
            return ibPointsPtr_();
        }

        List<vector>& getIbNormals()
        {
            return ibNormalsPtr_();
        }

        virtual void setIntpInfo() = 0;

        const DynamicLabelList& getSurfCells() const
        {
            return geom_model_->getSurfaceCellList()[Pstream::myProcNo()];
        }

        virtual void clearIntpInfo()
        {
            ibPointsPtr_.clear();
            ibNormalsPtr_.clear();
        }
    };
}
#endif
