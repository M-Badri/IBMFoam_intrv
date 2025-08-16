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
    nonConvex_body.C
\*---------------------------------------------------------------------------*/

#ifndef nonConvex_body_H
#define nonConvex_body_H

#include "stl_based.H"

namespace Foam
{
    class nonConvex_body
    : public stl_based
    {

        private:

        typedef DynamicList<label> DynamicLabelList;
        typedef DynamicList<scalar> DynamicScalarList;
        typedef DynamicList<vector> DynamicVectorList;
        typedef DynamicList<point> DynamicPointList;

        protected:


        public:

        // Constructors ------------------------------------------------
        nonConvex_body
        (
            const  fvMesh&   mesh,
            word      stlPath,
            scalar  thrSurf = 1
        )
        :
        stl_based
        (
            mesh,
            nonConvex,
            stlPath,
            thrSurf
        ){}

        // Destructors -------------------------------------------------
        ~nonConvex_body(){}

        void create_immersed_body
        (
            volScalarField& body,
            Field<label>& octreeField,
            List<labelList>& cellPoints
        ) final override;

        labelList getBBoxCellsByOctTree
        (
            label cellToCheck,
            bool& insideBB,
            vector& bBoxMin,
            vector& bBoxMax,
            List<DynamicLabelList>& bBoxCells,
            Field<label>& octreeField
        );

        std::shared_ptr<geom_model> getCopy() final override
        {
            std::shared_ptr<nonConvex_body> convex_bodyModel =
                std::make_shared<nonConvex_body>(mesh_,stlPath_,thrSurf_);

            convex_bodyModel->initVars(bodySurfMesh_.points(), thrSurf_, intSpan_, sdBasedLambda_);
            convex_bodyModel->setRhoS(rhoS_);
            convex_bodyModel->setM0(M0_);
            convex_bodyModel->getSurfaceCellList() = surfCells_;
            convex_bodyModel->getInternalCellList() = intCells_;

            return convex_bodyModel;
        }
    };
}

#endif
