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
    periodicBody.C
\*---------------------------------------------------------------------------*/

#ifndef periodicBody_H
#define periodicBody_H

#include "cluster_body.H"
#include "OFstream.H"

namespace Foam
{
    class periodicBody
    : public cluster_body
    {

        private:

        typedef DynamicList<label> DynamicLabelList;
        typedef DynamicList<scalar> DynamicScalarList;
        typedef DynamicList<vector> DynamicVectorList;
        typedef DynamicList<point> DynamicPointList;

        protected:


        public:

        // Constructors ------------------------------------------------
        periodicBody
        (
            const  fvMesh&   mesh,
            scalar  thrSurf = 1
        )
        :
        cluster_body
        (
            mesh,
            thrSurf
        )
        {}

        // Destructors -------------------------------------------------
        ~periodicBody(){}

 
        std::shared_ptr<geom_model> getCopy() final override
        {
            std::shared_ptr<periodicBody> periodicGeomModel = std::make_shared<periodicBody>
            (
                mesh_,
                thrSurf_
            );

            periodicGeomModel->ibGeomModelList = ibGeomModelList;
            return periodicGeomModel;
        }

        void getReferencedLists(
            List<DynamicLabelList>& intLists,
            List<DynamicLabelList>& surfLists,
            DynamicVectorList& referenceCoM
        ) final override;

        void bodyMovePoints(vector translVec) final override;

        void bodyScalePoints(scalar scaleFac) final override;

        void bodyRotatePoints
        (
            scalar rotAngle,
            vector axisOfRot
        ) final override;

        vector getCoM() final override;

        scalar getDC() final override;

       
        void setOwner() final override;

       
        label getOwner() final override;

     
        scalar& getM0() final override;

        vector getLVec(const point& toPoint) final override;

        void add_body_to_cluster(std::shared_ptr<geom_model> bodyToAdd)
        {
            ibGeomModelList.push_back(std::move(bodyToAdd));
        }

        bool shouldBeUnclustered();

        std::shared_ptr<geom_model> getRemGeomModel();

        void writeBodyInfo(label body_id, IOdictionary& dict, word path) final override
        {
            ibGeomModelList[0]->writeBodyInfo(body_id, dict, path);
        }
    };
}

#endif
