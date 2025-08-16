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
    stl_based.C

\*---------------------------------------------------------------------------*/

#ifndef stl_based_H
#define stl_based_H

#include "geom_model.H"

namespace Foam
{
    class stl_based
    : public geom_model
    {

        private:

        typedef DynamicList<label> DynamicLabelList;
        typedef DynamicList<scalar> DynamicScalarList;
        typedef DynamicList<vector> DynamicVectorList;
        typedef DynamicList<point> DynamicPointList;

        protected:
        triSurfaceMesh      bodySurfMesh_;
        word                stlPath_;
        autoPtr<triSurface> triSurf_;
        autoPtr<triSurfaceSearch> triSurfSearch_;

        pointField historyPoints_;

        void getIntersectionPoints
        (
            const label index,
            const treeBoundBox& cubeBb,
            DynamicPointList& intersectionPoints
        );

        void intersectBb
        (
            const point& p0,
            const point& p1,
            const point& p2,
            const treeBoundBox& cubeBb,
            DynamicPointList& intersectionPoints
        );

        public:

        // Constructors ------------------------------------------------
        stl_based
        (
            const  fvMesh&   mesh,
            const contactType cType,
            word      stlPath,
            scalar  thrSurf
        );

        // Destructors -------------------------------------------------
        virtual ~stl_based(){}

        virtual void create_immersed_body
        (
            volScalarField& body,
            Field<label>& octreeField,
            List<labelList>& cellPoints
        )
        {
            InfoH << iB_Info << "-- geom_modelMessage-- "
                << "notIplemented"<<endl;
        }

        virtual std::shared_ptr<geom_model> getCopy() = 0;

        vector add_modelReturnRandomPosition(
            const bool allActiveCellsInMesh,
            const boundBox  cellZoneBounds,
            Random&         randGen
        ) final override;

        void bodyMovePoints(vector translVec) final override;

        void bodyScalePoints(scalar scaleFac) final override;

        void bodyRotatePoints
        (
            scalar rotAngle,
            vector axisOfRot
        ) final override;

        void resetBody() final override
        {
            triSurfaceMesh tmpSurf(IOobject
                (
                    stlPath_,
                    mesh_,
                    IOobject::MUST_READ,
                    IOobject::NO_WRITE
                )
            );

            bodySurfMesh_.movePoints(tmpSurf.points());
        }

        vector getCoM() final override
        {
            pointField bSMeshPts = bodySurfMesh_.points();
            // get its center of mass
            vector CoM(vector::zero);
            forAll(bSMeshPts,point)
            {
                CoM += bSMeshPts[point];
            }
            CoM/= bSMeshPts.size();
            return CoM;
        }

        boundBox getBounds() final override
        {
            boundBox bounds(bodySurfMesh_.points(), false);
            return bounds;
        }

        void synchronPos(label owner = -1) final override;

        void writeBodyInfo(label body_id, IOdictionary& dict, word path) final override
        {
            triSurface triToRet(bodySurfMesh_);
            triToRet.write(path+"/stlFiles/"+name(body_id) + ".stl",".stl");
        }

        boolList pointInside(pointField pointI) final override
        {
            boolList returnList = triSurfSearch_().calcInside( pointI);
            return returnList;
        }

        bool pointInside(point pointI) final override
        {
            pointField points(1,pointI);
            boolList returnList = triSurfSearch_().calcInside( points);
            return returnList[0];
        }

        volumeType getVolumeType(sub_volume& sv, bool cIb)  final override;

        void getClosestPointAndNormal
        (
            const point& startPoint,
            const vector& span,
            point& closestPoint,
            vector& normal
        ) final override
        {
            pointIndexHit ibPointIndexHit = triSurfSearch_().nearest(startPoint, span);
            List<pointIndexHit> ibPointIndexHitList(1,ibPointIndexHit);
            vectorField normalVectorField;

            const triSurfaceMesh& ibTempMesh( bodySurfMesh_);
            ibTempMesh.getNormal(ibPointIndexHitList,normalVectorField);

            if(ibPointIndexHit.hit())
            {
                normal = normalVectorField[0];
                closestPoint = ibPointIndexHit.hitPoint();
            }
            else
            {
                InfoH << basic_Info << "Missing the closest point!" << endl;
                normal = startPoint - getCoM();
                closestPoint = getCoM();
            }
        }

        bool limitFinalSubVolume
        (
            const sub_volume& sv,
            bool cIb,
            boundBox& limBBox
        ) final override;

        void initVars(pointField surfPoints, scalar thrSurf, scalar intSpan, bool sdBasedLambda)
        {
            bodySurfMesh_.movePoints(surfPoints);
            triSurf_.reset(new triSurface(bodySurfMesh_));
            triSurfSearch_.reset(new triSurfaceSearch(triSurf_()));
            historyPoints_ = surfPoints;
            thrSurf_ = thrSurf;
            intSpan_ = intSpan;
            sdBasedLambda_ = sdBasedLambda;
        }

        void setCoM() final override
        {
            pointField bSMeshPts = bodySurfMesh_.points();
            CoM_    = vector::zero;
            forAll(bSMeshPts,point)
            {
                CoM_ += bSMeshPts[point];
            }
            CoM_/= bSMeshPts.size();
        }

        pointField getBodyPoints() final override
        {
            return bodySurfMesh_.points();
        }

        void setBodyPosition(pointField pos) final override;
    };
}
#endif
