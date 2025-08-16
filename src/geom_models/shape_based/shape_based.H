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
    shape_based.C

\*---------------------------------------------------------------------------*/

#ifndef shape_based_H
#define shape_based_H

#include "geom_model.H"
#include "meshSearch.H"

namespace Foam
{
    class shape_based
    : public geom_model
    {

        private:

        typedef DynamicList<label> DynamicLabelList;
        typedef DynamicList<scalar> DynamicScalarList;
        typedef DynamicList<vector> DynamicVectorList;
        typedef DynamicList<point> DynamicPointList;

        protected:

        public:

        // Constructors ------------------------------------------------
        shape_based
        (
            const  fvMesh&   mesh,
            const contactType cType,
            scalar  thrSurf
        );

        // Destructors -------------------------------------------------
        virtual ~shape_based(){}

        virtual void create_immersed_body
        (
            volScalarField& body,
            Field<label>& octreeField,
            List<labelList>& cellPoints
        )
        {
            InfoH << iB_Info << "-- shapeModelMessage-- "
                << "notIplemented"<<endl;
        }

        virtual std::shared_ptr<geom_model> getCopy() = 0;

        vector add_modelReturnRandomPosition(
            const bool allActiveCellsInMesh,
            const boundBox  cellZoneBounds,
            Random&         randGen
        ) final override;

        virtual void bodyMovePoints(vector translVec)
        {
            InfoH << iB_Info << "-- shapeModelMessage-- "
                << "notIplemented"<<endl;
        }

        virtual void bodyScalePoints(scalar scaleFac)
        {
            InfoH << iB_Info << "-- shapeModelMessage-- "
                << "notIplemented"<<endl;
        }

        virtual void bodyRotatePoints
        (
            scalar rotAngle,
            vector axisOfRot
        )
        {
            InfoH << iB_Info << "-- shapeModelMessage-- "
                << "notIplemented"<<endl;
        }

        virtual void resetBody()
        {
            InfoH << iB_Info << "-- shapeModelMessage-- "
                << "notIplemented"<<endl;
        }

        virtual vector getCoM()
        {
            InfoH << iB_Info << "-- shapeModelMessage-- "
                << "notIplemented"<<endl;
            return vector::zero;
        }

        virtual boundBox getBounds()
        {
            InfoH << iB_Info << "-- shapeModelMessage-- "
                << "notIplemented"<<endl;
            boundBox bound;
            return bound;
        }

        virtual void synchronPos(label owner)
        {
            InfoH << iB_Info << "-- shapeModelMessage-- "
                << "notIplemented"<<endl;
        }

        virtual void writeBodyInfo(label body_id, IOdictionary& dict, word path)
        {
            InfoH << iB_Info << "-- shapeModelMessage-- "
                << "notIplemented"<<endl;
        }

        virtual boolList pointInside(pointField pointI)
        {
            InfoH << iB_Info << "-- shapeModelMessage-- "
                << "notIplemented"<<endl;
            return {false};
        }


        virtual bool pointInside(point pointI)
        {
            InfoH << iB_Info << "-- geom_modelMessage-- "
                << "notIplemented"<<endl;
            return false;
        }

        virtual void getClosestPointAndNormal
        (
            const point& startPoint,
            const vector& span,
            point& closestPoint,
            vector& normal
        )
        {
            InfoH << iB_Info << "-- geom_modelMessage-- "
                << "notIplemented"<<endl;
        }

        virtual pointField sampleSurfacePoints()
        {
            InfoH << iB_Info << "-- shapeModelMessage-- "
                << "notIplemented"<<endl;
            pointField returnField;
            return returnField;
        }

        virtual void setCoM()
        {
            InfoH << iB_Info << "-- geom_modelMessage-- "
                << "notIplemented"<<endl;
        }
    };
}
#endif
