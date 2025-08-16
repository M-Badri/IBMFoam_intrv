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
    LS_interpol.C
\*---------------------------------------------------------------------------*/

#ifndef LS_interpol_H
#define LS_interpol_H

#include "dictionary.H"
#include "fvCFD.H"
#include "triSurfaceMesh.H"
#include "ib_interpol.H"
#include "LS_interpolInfo.H"

namespace Foam
{
    class LS_interpol
    : public ib_interpol
    {
        private:
        typedef DynamicList<label> DynamicLabelList;
        typedef DynamicList<scalar> DynamicScalarList;
        typedef DynamicList<vector> DynamicVectorList;
        typedef DynamicList<point> DynamicPointList;

        const  fvMesh&   mesh_;

        // values corresponding to cell cell searching
        const  scalar   distFactor_;
        const  scalar   radiusFactor_;
        const  scalar   angleFactor_;
        const  scalar   maxCCRows_;

        void imposeDirichletCondition
        (
            LS_interpolInfo& intpInfo,
            volVectorField& Ui,
            vectorField& ibPointsVal,
            const Foam::fvMesh& mesh
        );

        void adjustPhi
        (
            LS_interpolInfo& intpInfo,
            surfaceScalarField& phi
        );

        public:

        // Constructors ------------------------------------------------
        LS_interpol(
            const fvMesh& mesh,
            scalar distFactor,
            scalar radiusFactor,
            scalar angleFactor,
            scalar maxCCRows
        );

        // Destructors -------------------------------------------------
        ~LS_interpol();

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
            return new LS_interpolInfo
            (
                mesh,
                gModel,
                distFactor_,
                radiusFactor_,
                angleFactor_,
                maxCCRows_
            );
        };
    };
}
#endif
