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
    class for clustering of imm_bodies

SourceFiles
    cluster_body.C
\*---------------------------------------------------------------------------*/

#ifndef cluster_body_H
#define cluster_body_H

#include "geom_model.H"
#include "OFstream.H"

namespace Foam
{
    class cluster_body
    : public geom_model
    {

        private:

        typedef DynamicList<label> DynamicLabelList;
        typedef DynamicList<scalar> DynamicScalarList;
        typedef DynamicList<vector> DynamicVectorList;
        typedef DynamicList<point> DynamicPointList;

        protected:

        std::vector<std::shared_ptr<geom_model>> ibGeomModelList;

        public:

        // Constructors ------------------------------------------------
        cluster_body
        (
            const  fvMesh&   mesh,
            scalar  thrSurf = 1
        )
        :
        geom_model
        (
            mesh,
            cluster,
            thrSurf
        )
        {}

        // Destructors -------------------------------------------------
        ~cluster_body(){}

        // Public subclasses -------------------------------------------

        // Public variables --------------------------------------------

        // Public methods ----------------------------------------------

        void create_immersed_body
        (
            volScalarField& body,
            Field<label>& octreeField,
            List<labelList>& cellPoints
        ) final override;

        virtual std::shared_ptr<geom_model> getCopy() = 0;

        // Connect int/surf cells with correct center of mass
        // for fluid-solid interaction
        virtual void getReferencedLists(
            List<DynamicLabelList>& intLists,
            List<DynamicLabelList>& surfLists,
            DynamicVectorList& referenceCoM
        ) = 0;

        void updateSurfList();

        virtual List<DynamicLabelList>& getSurfaceCellList() final override
        {
            updateSurfList();
            return surfCells_;
        };

        void updateIntList();

        virtual List<DynamicLabelList>& getInternalCellList() final override
        {
            updateIntList();
            return intCells_;
        };

        void calculateGeometricalProperties(
            volScalarField& body
        ) final override;

        void calculateGeometricalPropertiesParallel(
            volScalarField& body
        ) final override;

        void setMassAndInertia();

        virtual void bodyMovePoints(vector translVec) = 0;

        virtual void bodyScalePoints(scalar scaleFac) = 0;

        virtual void bodyRotatePoints
        (
            scalar rotAngle,
            vector axisOfRot
        ) = 0;

        virtual vector getCoM() = 0;

        boundBox getBounds() final override;

        void synchronPos(label owner = -1) final override;

        boolList pointInside(pointField pointI) final override;

        bool pointInside(point pointI) final override;

        // Get characteristic particle diameter
        virtual scalar getDC() = 0;

        void getClosestPointAndNormal
        (
            const point& startPoint,
            const vector& span,
            point& closestPoint,
            vector& normal
        ) final override;


        virtual label getOwner() = 0;

    
        virtual scalar& getM0() = 0;

        virtual vector getLVec(const point& toPoint) = 0;

        void resetBody(volScalarField& body) final override;

        bool isCluster() const final override
        {
            return true;
        }

        std::vector<std::shared_ptr<geom_model>>& get_cluster_bodies()
        {
            return ibGeomModelList;
        }

        List<std::shared_ptr<boundBox>> getBBoxes() final override;

        virtual void writeBodyInfo(label body_id, IOdictionary& dict, word path)
        {
            InfoH << iB_Info << "-- geom_modelMessage-- "
                << "custer body writeBodyInfo notIplemented"<<endl;
        }
    };
}

#endif
