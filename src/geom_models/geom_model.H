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
    geom_model.C
\*---------------------------------------------------------------------------*/

#ifndef geom_model_H
#define geom_model_H

#include <vector>
#include <memory>

#include "dictionary.H"
#include "fvCFD.H"
#include "triSurfaceMesh.H"

#include "contactType.H"
#include "sub_volume.H"
#include "verletPoint.H"

#include "outputIBMFoam.H"
using namespace outputIBMFoam;

namespace Foam
{

    template<class Type>
    Foam::tmp<Foam::Field<Type>>
    filterField
    (
        const Field<Type>& field,
        const labelList& cellIDs
    )
    {
        if (isNull(cellIDs))
        {
            return field;
        }
        else
        {
            return tmp<Field<Type>>(new Field<Type>(field, cellIDs));
        }
    }

    class geom_model
    {

        private:

        typedef DynamicList<label> DynamicLabelList;
        typedef DynamicList<scalar> DynamicScalarList;
        typedef DynamicList<vector> DynamicVectorList;
        typedef DynamicList<point> DynamicPointList;

        const contactType contactType_;

        protected:

        const  fvMesh&   mesh_;


        List<DynamicLabelList>  surfCells_;
        List<DynamicLabelList>  intCells_;
        DynamicLabelList ibPartialVolume_;
        label owner_;

        label cellToStartInCreateIB_;
        scalar  thrSurf_;

        scalar  intSpan_;
        bool    sdBasedLambda_;
        boundBox curMeshBounds_;
        autoPtr<HashTable<const labelList&, label, Hash<label>>> cachedNeighbours_;

        scalar M_;
        scalar M0_;                                                  
        label nCells_;


        vector CoM_;

     
        symmTensor  I_;

        std::shared_ptr<boundBox> bBox_;

    
        scalar  dC_;

      
        dimensionedScalar rhoS_;

        public:

        // Constructors ------------------------------------------------
        geom_model
        (
            const  fvMesh&   mesh,
            const contactType cType,
            scalar  thrSurf = 1
        );

        // Destructors -------------------------------------------------
        virtual ~geom_model();

        virtual void create_immersed_body
        (
            volScalarField& body,
            Field<label>& octreeField,
            List<labelList>& cellPoints
        )
        {
            InfoH << iB_Info << "1-- geom_modelMessage-- "
                << "notIplemented"<<endl;
        }

        virtual std::shared_ptr<geom_model> getCopy() = 0;

        virtual List<DynamicLabelList>& getSurfaceCellList()
        {
            return surfCells_;
        };

        virtual void getReferencedLists(
            List<DynamicLabelList>& intLists,
            List<DynamicLabelList>& surfLists,
            DynamicVectorList& referenceCoM
        )
        {
            intLists.resize(1, intCells_[Pstream::myProcNo()]);
            surfLists.resize(1, surfCells_[Pstream::myProcNo()]);
            referenceCoM.resize(1, getCoM());
        };

        virtual List<DynamicLabelList>& getInternalCellList()
        {
            return intCells_;
        };

        virtual vector add_modelReturnRandomPosition(
            const bool allActiveCellsInMesh,
            const boundBox  cellZoneBounds,
            Random&         randGen
        )
        {
            InfoH << iB_Info << "2-- geom_modelMessage-- "
                << "notIplemented"<<endl;
            return vector::zero;
        }

        virtual void bodyMovePoints(vector translVec)
        {
            InfoH << iB_Info << "3-- geom_modelMessage-- "
                << "notIplemented"<<endl;
        }

        virtual void bodyScalePoints(scalar scaleFac)
        {
            InfoH << iB_Info << "4-- geom_modelMessage-- "
                << "notIplemented"<<endl;
        }

        virtual void bodyRotatePoints
        (
            scalar rotAngle,
            vector axisOfRot
        )
        {
            InfoH << iB_Info << "5-- geom_modelMessage-- "
                << "notIplemented"<<endl;
        }

        virtual void resetBody()
        {
            InfoH << iB_Info << "6-- geom_modelMessage-- "
                << "notIplemented"<<endl;
        }

        virtual vector getCoM()
        {
            InfoH << iB_Info << "7-- geom_modelMessage-- "
                << "notIplemented"<<endl;
            return vector::zero;
        }

        virtual boundBox getBounds()
        {
            InfoH << iB_Info << "8-- geom_modelMessage-- "
                << "notIplemented"<<endl;
            boundBox bound;
            return bound;
        }

        virtual void synchronPos(label owner = -1)
        {
            InfoH << iB_Info << "9-- geom_modelMessage-- "
                << "notIplemented"<<endl;
        }

        virtual void writeBodyInfo(label body_id, IOdictionary& dict, word path)
        {
            InfoH << iB_Info << "10-- geom_modelMessage-- "
                << "notIplemented"<<endl;
        }

        virtual boolList pointInside(pointField pointI)
        {
            InfoH << iB_Info << "11-- geom_modelMessage-- "
                << "notIplemented"<<endl;
            return {false};
        }

        virtual bool pointInside(point pointI)
        {
            InfoH << iB_Info << "12-- geom_modelMessage-- "
                << "notIplemented"<<endl;
            return false;
        }

        virtual volumeType getVolumeType(sub_volume& sv, bool cIb)
        {
            InfoH << iB_Info << "13-- geom_modelMessage-- "
                << "notIplemented"<<endl;
            return volumeType::unknown;
        }

        virtual scalar getDC()
        {
            return dC_;
        }

        virtual void getClosestPointAndNormal
        (
            const point& startPoint,
            const vector& span,
            point& closestPoint,
            vector& normal
        )
        {
            InfoH << iB_Info << "14-- geom_modelMessage-- "
                << "notIplemented"<<endl;
        }

        virtual bool limitFinalSubVolume
        (
            const sub_volume& sv,
            bool cIb,
            boundBox& limBBox
        )
        {
            InfoH << iB_Info << "15-- geom_modelMessage-- "
                << "notIplemented"<<endl;
            return false;
        }

        virtual void setCoM()
        {
            InfoH << iB_Info << "16-- geom_modelMessage-- "
                << "notIplemented"<<endl;
        }

        constexpr contactType getcType()
        {return contactType_;}

        virtual void setOwner()
        {
            Pstream::gatherList(ibPartialVolume_, 0);
            Pstream::scatter(ibPartialVolume_, 0);
            for (label i = 0; i < ibPartialVolume_.size(); i++)
            {
                if (ibPartialVolume_[i] == max(ibPartialVolume_))
                {
                    //set owner of the IB which will move this IB
                    owner_ = i;
                    break;
                }
            }
        }

        virtual label getOwner()
        {
            return owner_;
        }

        void setIntSpan(scalar intSpan){intSpan_ = intSpan;}

        void setSdBasedLambda(bool sdBasedLambda){sdBasedLambda_ = sdBasedLambda;}

        void resetHashTable(){cachedNeighbours_.clear();}

        virtual void setRhoS(dimensionedScalar rhoS){rhoS_ = rhoS;}
        dimensionedScalar getRhoS(){return rhoS_;}

        virtual void calculateGeometricalProperties(
            volScalarField& body
        );

        virtual void calculateGeometricalPropertiesParallel(
            volScalarField& body
        );

        void addToMAndI
        (
            volScalarField& body,
            DynamicLabelList& labelCellLst
        );

        void compute_body_charPars();

        virtual scalar& getM0()
        {
            return M0_;
        }

        void setM0(scalar M0)
        {
            M0_ = M0;
        }

        scalar& getM()
        {
            return M_;
        }

        void setM(scalar M)
        {
            M_ = M;
        }


        label& getNCells()
        {
            return nCells_;
        }

        void setNCells(label nCells)
        {
            nCells_ = nCells;
        }

        const dimensionedScalar getRhoS()    const
        {
            return rhoS_;
        }

        virtual vector getLVec(const point& toPoint)
        {
            return toPoint - getCoM();
        }

        const symmTensor& getI()    const
        {
            return I_;
        }

        void setI(const symmTensor& I)
        {
            I_ = I;
        }

        virtual void resetBody(volScalarField& body);

        virtual bool isCluster() const
        {
            return false;
        }

        virtual bool checkCluster()
        {
            InfoH << iB_Info << "17-- geom_modelMessage-- "
                << "not cluster model"<<endl;
                return false;
        }

        bool isBBoxInMesh();

        DynamicLabelList getPotentSurfCells
        (
            volScalarField& body,
            HashTable<bool, label, Hash<label>>& cellInside,
            List<labelList>& cellPoints
        );

        void correctSurfCells
        (
            volScalarField& body,
            DynamicLabelList& potentSurfCells,
            HashTable<bool, label, Hash<label>>& cellInside,
            List<labelList>& cellPoints
        );

        virtual List<std::shared_ptr<boundBox>> getBBoxes();

        virtual pointField getBodyPoints()
        {
            InfoH << iB_Info << "18-- geom_modelMessage-- "
                << "not cluster model"<<endl;
            pointField emptyField;
            return emptyField;
        }

        virtual  void setBodyPosition(pointField pos)
        {
            InfoH << iB_Info << "22-- geom_modelMessage-- "
                << "notIplemented"<<endl;
        }

    };
}
#endif
