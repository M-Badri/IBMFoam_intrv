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
    sphere_body.C

\*---------------------------------------------------------------------------*/

#ifndef sphere_body_H
#define sphere_body_H

#include "shape_based.H"
#include "OFstream.H"

#include <assert.h>

namespace Foam
{
    class sphere_body
    : public shape_based
    {

        private:

        typedef DynamicList<label> DynamicLabelList;
        typedef DynamicList<scalar> DynamicScalarList;
        typedef DynamicList<vector> DynamicVectorList;
        typedef DynamicList<point> DynamicPointList;

        vector position_;
        scalar radius_;
        scalar sphereM_;

        vector startPosition_;
        scalar startRadius_;

        vector historyPosition_;
        scalar historyRadius_;

        public:

        // Constructors ------------------------------------------------
        sphere_body
        (
            const  fvMesh&   mesh,
            vector startPos,
            scalar radius,
            scalar  thrSurf = 1
        )
        :
        shape_based
        (
            mesh,
            sphere,
            thrSurf
        ),
        position_(startPos),
        radius_(radius),
        sphereM_((4*Foam::constant::mathematical::pi*pow(radius,3))/3)
        {
            startPosition_ = position_;
            startRadius_ = radius_;

            historyPosition_ = position_;
            historyRadius_ = radius_;
        }

        // Destructors -------------------------------------------------
        ~sphere_body(){}


        void create_immersed_body
        (
            volScalarField& body,
            Field<label>& octreeField,
            List<labelList>& cellPoints
        ) final override;

        label getCellInBody
        (
            Field<label>& octreeField
        );

        std::shared_ptr<geom_model> getCopy() final override
        {
            std::shared_ptr<sphere_body> sphere_bodyModel =
                std::make_shared<sphere_body>(mesh_,position_,radius_,thrSurf_);

            sphere_bodyModel->initVars(intSpan_,sdBasedLambda_);
            sphere_bodyModel->setRhoS(rhoS_);
            sphere_bodyModel->setM0(M0_);
            sphere_bodyModel->getSurfaceCellList() = surfCells_;
            sphere_bodyModel->getInternalCellList() = intCells_;

            return sphere_bodyModel;
        }

        void bodyMovePoints(vector translVec) final override
        {position_ += translVec;}

        void bodyScalePoints(scalar scaleFac) final override
        {radius_ *= scaleFac;}

        void bodyRotatePoints
        (
            scalar rotAngle,
            vector axisOfRot
        ) final override
        {}

        void resetBody() final override
        {
            position_ = startPosition_;
            radius_ = startRadius_;
        }

        vector getCoM() final override
        {return position_;}

        boundBox getBounds() final override
        {
            vector one(1,1,1);
            one *= radius_;
            boundBox bounds(position_-one,position_+one);
            return bounds;
        }

        void synchronPos(label owner = -1) final override;

        void writeBodyInfo(label body_id, IOdictionary& dict, word path) final override
        {
            dictionary sphereDict("sphere");
            sphereDict.add("position", position_);
            sphereDict.add("radius", radius_);

            dict.add("sphere",sphereDict);
        }

        boolList pointInside(pointField pointI) final override;


        bool pointInside(point pointI) final override;

        void getClosestPointAndNormal
        (
            const point& startPoint,
            const vector& span,
            point& closestPoint,
            vector& normal
        ) final override
        {
            if(startPoint == position_)
            {
                closestPoint = startPoint;
                normal = vector::zero;
            }
            else
            {
                normal = startPoint-position_;
                normal /= mag(normal);
                closestPoint = position_ + normal * radius_;
            }
        }

        void initVars(scalar intSpan, bool sdBasedLambda)
        {
            intSpan_ = intSpan;
            sdBasedLambda_ = sdBasedLambda;
        }

        pointField sampleSurfacePoints() final override;

        scalar getDC() final override
        {
            return 2*radius_;
        }

        scalar& getM0() final override
        {
            if(case3D)
            {
                return sphereM_;
            }
            else
            {
                return M0_;
            }
        }

        void setRhoS(dimensionedScalar rhoS)
        {
            rhoS_ = rhoS;
            sphereM_ = ((4*Foam::constant::mathematical::pi*pow(radius_,3))/3)*rhoS_.value();
        }

        void setCoM() final override
        {
            CoM_ = position_;
        }

        void setBodyPosition(pointField pos) final override
        {
            assert(pos.size() == 1 && "Only one point is allowed for sphere body");
            position_ = pos[0];
        }

        pointField getBodyPoints() final override
        {
            pointField points(1,position_);
            return points;
        }
    };
}

#endif
