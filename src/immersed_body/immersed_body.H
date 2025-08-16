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
\*---------------------------------------------------------------------------*/

#ifndef immersed_body_H
#define immersed_body_H

#include "dictionary.H"
#include "coordinateSystem.H"
#include "coordinateSystems.H"
#include "wordList.H"
#include "labelList.H"
#include "dimensionedScalar.H"
#include "dimensionedTensor.H"
#include "primitiveFieldsFwd.H"
#include "volFieldsFwd.H"
#include "fvMatricesFwd.H"
#include "triSurfaceMesh.H"

#include "geom_model.H"
#include "stl_based.H"
#include "convex_body.H"
#include "nonConvex_body.H"

#include "ibContactClass.H"
#include "wallContactInfo.H"
#include "interpol_info.H"
#include "ib_interpol.H"

#include "virtual_mesh.H"
namespace Foam
{

    class immersed_body
    {

        private:

        // Private data types ------------------------------------------

        typedef DynamicList<label> DynamicLabelList;
        typedef DynamicList<scalar> DynamicScalarList;
        typedef DynamicList<vector> DynamicVectorList;
        typedef DynamicList<point> DynamicPointList;

        word    body_name_;
        bool    isActive_;
        dictionary&             immersedDict_;
        const  fvMesh&   mesh_;
        dictionary&             trans_properties_;

        std::shared_ptr<geom_model> geom_model_;

        autoPtr<interpol_info> intpInfo_;

        autoPtr<ibContactClass> ibContactClass_;

        autoPtr<ibContactVars> ibCntVars_;

        autoPtr<wallContactInfo> wallCntInfo_;

        List<labelList>& cellPoints_;

        vector Axis_;
        vector AxisOld_;

        scalar omega_;
        scalar omegaOld_;

        vector Vel_;
        vector VelOld_;
        vector a_;

        vector alpha_;

        vector totalAngle_;

        forces FCoupling_;

        forces FContact_;

        scalar CoNum_;

        scalar meanCoNum_;

        dimensionedScalar rhoF_;
        label   body_id_;

        bool    update_torq_;

        int     bodyOperation_;

        void  updateImmersedBody
        (
            volScalarField& body,
            volVectorField& f
        );

        void updateCoupling
        (
            volScalarField& body,
            volVectorField& f
        );

        Field<label> octreeField_;
        label cellToStartInCreateIB_;

        bool    startSynced_;

        tensor  totRotMatrix_;
        bool    sdBasedLambda_;
        scalar  intSpan_;
        scalar  charCellSize_;                                          
        label   refineBuffers_;                                         
        label   recomputeM0_;                                          
        label   t_to_set_static_;
        vector  staticContactPost_;
        protected:


        public:

        immersed_body
        (
            word body_name,
            const Foam::fvMesh& mesh,
            dictionary& IBMFoamDict,
            dictionary& transportProperties,
            label body_id,
            label recomputeM0,
            std::shared_ptr<geom_model> b_geomModel,
            autoPtr<ib_interpol>& ibIntp,
            List<labelList>& cellPoints
        );

        // Destructors -------------------------------------------------
        ~immersed_body();

        // Public methods ----------------------------------------------

        // move immersed body
        void moveImmersedBody
        (
            scalar deltaT = -1.0
        );
        
        void printBodyInfo();
        
        void updateMovement
        (
            scalar deltaT
        );
        void updateMovement
        (
            vector Vel,
            vector Axis,
            scalar omega
        );

        void updateMovementComp
        (
            scalar deltaT,
            vector Vel,
            vector Axis,
            scalar omega
        );

        void create_immersed_body
        (
            volScalarField& body,
            volScalarField& refine_F,
            bool synchCreation = true
        );

        void syncCreateImmersedBody
        (
            volScalarField& body,
            volScalarField& refine_F
        );

        void syncImmersedBodyParralell1
        (
            volScalarField& body,
            volScalarField& refine_F
        );

        void syncImmersedBodyParralell2
        (
            volScalarField& body,
            volScalarField& refine_F
        );

        void constructRefineField
        (
            volScalarField& body,
            volScalarField& refine_F,
            DynamicLabelList cellsToIterate,
            DynamicLabelList startLevel
        );

        // getters

        interpol_info& getIntpInfo()
        {
            return intpInfo_();
        }

        void clearIntpInfo()
        {
            if(intpInfo_.valid())
            {
                intpInfo_().clearIntpInfo();
            }
        }

        const bool& checkWallContact() const
        {
            return ibContactClass_->checkWallContact();
        }

        const label& getBodyId() const
        {
            return body_id_;
        }

        const vector& getVel()  const
        {
            return Vel_;
        }

        const scalar& getOmega()    const
        {
            return omega_;
        }

        const vector& getAxis()     const
        {
            return Axis_;
        }

        const scalar& getCoNum()   const
        {
            return CoNum_;
        }

        const scalar& getMeanCoNum()   const
        {
            return meanCoNum_;
        }

        scalar getDC()
        {
            return geom_model_->getDC();
        }

        const int& getbodyOperation()   const
        {
            return bodyOperation_;
        }

        List<DynamicLabelList>& getSurfaceCellList()
        {
            return geom_model_->getSurfaceCellList();
        };

        List<DynamicLabelList>& getInternalCellList()
        {
            return geom_model_->getInternalCellList();
        };

        void postContactUpdateBodyField
        (
            volScalarField& body,
            volScalarField& refine_F
        );

        void recreateBodyField(volScalarField& body, volScalarField& refine_F);

        void update_vector_field
        (
            volVectorField& VS,
            word VName,
            volScalarField& body
        );

        void updateContactForces
        (
            const forces& F
        )
        {
            FContact_ += F;
        }

        void setWallContact(bool inContact)
        {
            ibContactClass_->setWallContact(inContact);
        }

        void post_pimple_update_immersed_body
        (
            volScalarField& body,
            volVectorField& f
        );

        void computeBodyCoNumber();

        void printMomentum();

        void printStats();

        bool& retUpdateTorque()
        {
            return update_torq_;
        }

        void reset_contact_forces()
        {
            FContact_.clear();
        }

        void switchActiveOff(volScalarField& body);
        const bool& get_is_active()
        {
            return isActive_;
        }

        const bool& get_start_synced() const
        {
            return startSynced_;
        }

        void initSyncWithFlow(const volVectorField& U);

        void pimpleUpdate
        (
            volScalarField& body,
            volVectorField& f
        );

        const bool& getSDBasedLambda() const
        {
            return sdBasedLambda_;
        }

        const scalar& getIntSpan() const
        {
            return intSpan_;
        }

        void compute_body_charPars(){geom_model_->compute_body_charPars();}

        void check_if_in_domain(volScalarField& body);

        const label& getrecomputeM0() const
        {
            return recomputeM0_;
        }
        void recomputedM0()
        {
            recomputeM0_--;
        }

        void inContact_with_static(bool inContact)
        {
            ibContactClass_->inContact_with_static(inContact);
        }

        void chceckBodyOp();

        void setRestartSim
        (
            vector vel,
            scalar angVel,
            vector axisRot,
            bool setStatic,
            label timesInContact
        );

        geom_model& get_geom_model()
        {
            return *geom_model_;
        }

        std::shared_ptr<geom_model>& get_geom_modelPtr()
        {
            return geom_model_;
        }

        ibContactClass& get_ib_contact_class()
        {
            return ibContactClass_();
        }

        ibContactVars& getContactVars()
        {
            return ibCntVars_();
        }

        wallContactInfo& get_wall_cnt_info()
        {
            return wallCntInfo_();
        }

        void record_body_info(IOdictionary& dict, word path)
        {
            dict.add("body_id",body_id_);
            dict.add("body_name",body_name_);
            dict.add("Vel",Vel_);
            dict.add("omega",omega_);
            dict.add("Axis",Axis_);
            dict.add("static",bodyOperation_ == 0);
            dict.add("timeStepsInContWStatic", ibContactClass_
                ->getTimeStepsInContWStatic());
            geom_model_->writeBodyInfo(body_id_,dict, path);
        }


        void update_old_movement_vars()
        {
            VelOld_   = Vel_;
            AxisOld_  = Axis_;
            omegaOld_ = omega_;
        }

        void writeIbInfo()
        {
            InfoH << iB_Info;
            InfoH << "-- body " << body_id_ << " center of mass       : " << geom_model_->getCoM() << endl;
            InfoH << "-- body " << body_id_ << " linear velocity      : " << Vel_ << endl;
            InfoH << "-- body " << body_id_ << " angluar velocity     : " << omega_ << endl;
            InfoH << "-- body " << body_id_ << " axis of rotation     : " << Axis_ << endl;
        }

        vectorField getUatIbPoints();
    };
}
#endif
