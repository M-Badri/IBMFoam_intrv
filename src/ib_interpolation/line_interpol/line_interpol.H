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
    class for addition and removal of imm_bodies

SourceFiles
    line_interpol.C
\*---------------------------------------------------------------------------*/

#ifndef line_interpol_H
#define line_interpol_H

#include "dictionary.H"
#include "fvCFD.H"
#include "triSurfaceMesh.H"
#include "interpolationCellPoint.H"
#include "interpolationCell.H"

#include "ib_interpol.H"
#include "line_interpolInfo.H"

namespace Foam
{
    class line_interpol
    : public ib_interpol
    {
        private:
        typedef DynamicList<label> DynamicLabelList;
        typedef DynamicList<scalar> DynamicScalarList;
        typedef DynamicList<vector> DynamicVectorList;
        typedef DynamicList<point> DynamicPointList;

        autoPtr<interpolation<vector>> interpV_;
        dictionary& interpDict_;

        void correctVelocity
        (
            line_interpolInfo& intpInfo,
            volVectorField& Ui,
            vectorField& ibPointsVal,
            const Foam::fvMesh& mesh
        );

        void getCurVelocity
        (
            List<List<interpol_point>>& interpol_points
        );

        List<label> getIntOrder
        (
            List<List<interpol_point>>& interpol_points
        );

        public:

        // Constructors ------------------------------------------------
        line_interpol(dictionary& interpDict);

        // Destructors -------------------------------------------------
        ~line_interpol();

        void resetInterpolator(volVectorField& V) override
        {
            interpV_ = interpolation<vector>::New(interpDict_, V);
        }

        void ib_interpolate
        (
            interpol_info& intpInfo,
            volVectorField& Ui,
            vectorField ibPointsVal,
            const Foam::fvMesh& mesh
        ) override;

        interpol_info* createIntpInfo
        (
            const Foam::fvMesh& mesh,
            std::shared_ptr<geom_model>& gModel
        ) override
        {
            return new line_interpolInfo
            (
                mesh,
                gModel
            );
        };
    };
}
#endif
