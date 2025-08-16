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
#include "wallSubContactInfo.H"

#include "inter_adhesion.H"
#include "wall_mat_info.H"

#include "virtual_meshLevel.H"
#include "wall_plane_info.H"
#include "contact_model_info.H"
using namespace Foam;
//---------------------------------------------------------------------------//
wallSubContactInfo::wallSubContactInfo
(
    List<Tuple2<point,boundBox>> contactBBData,
    List<Tuple2<point,boundBox>> planeBBData,
    List<string> contactPatches,
    List<Tuple2<point,boundBox>> internalBBData,
    HashTable<physicalProperties,string,Hash<string>> wallMeanPars,
    boundBox BB,
    label body_id
)
:
contactPatches_(contactPatches),
internalBBData_(internalBBData),
wallMeanPars_(wallMeanPars),
BB_(BB),
body_id_(body_id)
{
    forAll(contactBBData,cBD)
    {
        vector sub_volumeNVector = vector(
            floor((contactBBData[cBD].second().span()[0]/virtual_meshLevel::getCharCellSize())*virtual_meshLevel::getLevelOfDivision()),
            floor((contactBBData[cBD].second().span()[1]/virtual_meshLevel::getCharCellSize())*virtual_meshLevel::getLevelOfDivision()),
            floor((contactBBData[cBD].second().span()[2]/virtual_meshLevel::getCharCellSize())*virtual_meshLevel::getLevelOfDivision())
        );
        if(cmptMin(sub_volumeNVector)<SMALL)
        {
     
            for(int i=0;i<3;i++)
            {
                if(sub_volumeNVector[i] <SMALL)
                {                
                    sub_volumeNVector[i] = 1;
                    contactBBData[cBD].second().min()[i] -=virtual_meshLevel::getCharCellSize()/virtual_meshLevel::getLevelOfDivision()*0.5;
                    contactBBData[cBD].second().max()[i] +=virtual_meshLevel::getCharCellSize()/virtual_meshLevel::getLevelOfDivision()*0.5;
                }  
            }
        
        }
     

        autoPtr<virtual_meshWallInfo> vmWInfo(
            new virtual_meshWallInfo(
                contactBBData[cBD].second(),
                contactBBData[cBD].first(),
                sub_volumeNVector,
                virtual_meshLevel::getCharCellSize(),
                pow(virtual_meshLevel::getCharCellSize()/virtual_meshLevel::getLevelOfDivision(),3)
            )  
        );
        vmWInfoList_.append(vmWInfo);
    }

    forAll(planeBBData,pBD)
    {
        vector sub_volumeNVector = vector(
            ceil((planeBBData[pBD].second().span()[0]/virtual_meshLevel::getCharCellSize()))*virtual_meshLevel::getLevelOfDivision(),
            ceil((planeBBData[pBD].second().span()[1]/virtual_meshLevel::getCharCellSize()))*virtual_meshLevel::getLevelOfDivision(),
            ceil((planeBBData[pBD].second().span()[2]/virtual_meshLevel::getCharCellSize()))*virtual_meshLevel::getLevelOfDivision()
        );

        for(int i=0;i<3;i++)
        {
            if(sub_volumeNVector[i] == planeBBData[pBD].second().minDim())
            {                
                sub_volumeNVector[i] = 1;
                planeBBData[pBD].second().min()[i] -=virtual_meshLevel::getCharCellSize()/virtual_meshLevel::getLevelOfDivision()*0.5;
                planeBBData[pBD].second().max()[i] +=virtual_meshLevel::getCharCellSize()/virtual_meshLevel::getLevelOfDivision()*0.5;
            }  
        }

        autoPtr<virtual_meshWallInfo> vmWInfo(
            new virtual_meshWallInfo(
                planeBBData[pBD].second(),
                planeBBData[pBD].first(),
                sub_volumeNVector,
                virtual_meshLevel::getCharCellSize(),
                pow(virtual_meshLevel::getCharCellSize()/virtual_meshLevel::getLevelOfDivision(),3)
            )
        );
        vmPlaneInfoList_.append(vmWInfo);
    }
}

wallSubContactInfo::~wallSubContactInfo()
{
}
//---------------------------------------------------------------------------//
vector wallSubContactInfo::getLVec(wallContactVars& wallCntvar, ibContactClass ibCClass)
{
    return ibCClass.get_geom_model().getLVec(wallCntvar.contactCenter_);
}
//---------------------------------------------------------------------------//
vector wallSubContactInfo::getVeli(wallContactVars& wallCntvar, ibContactVars& cVars)
{
    return (-((wallCntvar.lVec_-cVars.Axis_
        *((wallCntvar.lVec_) & cVars.Axis_))
        ^ cVars.Axis_)*cVars.omega_+ cVars.Vel_);
}
//---------------------------------------------------------------------------//
void wallSubContactInfo::evalVariables(
    wallContactVars& wallCntvar,
    ibContactClass& ibCClass,
    ibContactVars& cVars
)
{
    reduceM_ =
    (
        ibCClass.get_geom_model().getM0()
        *ibCClass.get_geom_model().getM0()
        /(ibCClass.get_geom_model().getM0()
        +ibCClass.get_geom_model().getM0())
    );

    wallCntvar.lVec_ = getLVec(wallCntvar,ibCClass);
    wallCntvar.Veli_ = getVeli(wallCntvar, cVars);

    wallCntvar.Vn_ = -(wallCntvar.Veli_ - vector::zero) & wallCntvar.contactNormal_;
    wallCntvar.Lc_ = (contact_model_info::getLcCoeff())*mag(wallCntvar.lVec_)*mag(wallCntvar.lVec_)/(mag(wallCntvar.lVec_) + mag(wallCntvar.lVec_));
    

    wallCntvar.curAdhN_ = min
    (
        wallCntvar.getMeanCntPar().maxAdhN_,
        max(wallCntvar.curAdhN_, wallCntvar.getMeanCntPar().aY_
            *wallCntvar.contactVolume_
            /(sqr(wallCntvar.Lc_)*8*Foam::constant::mathematical::pi))
    );
}
//---------------------------------------------------------------------------//
vector wallSubContactInfo::getFNe(wallContactVars& wallCntvar)
{
    return (wallCntvar.getMeanCntPar().aY_*wallCntvar.contactVolume_
        /(wallCntvar.Lc_+SMALL))*wallCntvar.contactNormal_;
}
//---------------------------------------------------------------------------//
vector wallSubContactInfo::getFA(wallContactVars& wallCntvar)
{
    return ((sqrt(8*Foam::constant::mathematical::pi
        *wallCntvar.getMeanCntPar().aY_
        *wallCntvar.curAdhN_*wallCntvar.contactVolume_))
        *wallCntvar.contactNormal_);
}
//---------------------------------------------------------------------------//
vector wallSubContactInfo::getFNd(wallContactVars& wallCntvar)
{
    physicalProperties& meanCntPar(wallCntvar.getMeanCntPar());
    return ((meanCntPar.reduceBeta_*sqrt(meanCntPar.aY_
            *reduceM_*wallCntvar.contactArea_/(wallCntvar.Lc_+SMALL))*
            wallCntvar.Vn_)*wallCntvar.contactNormal_);

}
//---------------------------------------------------------------------------//
vector wallSubContactInfo::getFt(wallContactVars& wallCntvar, scalar deltaT)
{
    physicalProperties& meanCntPar(wallCntvar.getMeanCntPar());
    // project last Ft into a new direction
    vector FtLastP(wallCntvar.FtPrev_
        - (wallCntvar.FtPrev_ & wallCntvar.contactNormal_)
        *wallCntvar.contactNormal_);
    // scale projected Ft to have same magnitude as FtLast
    vector FtLastS(mag(wallCntvar.FtPrev_) * (FtLastP/(mag(FtLastP)+SMALL)));
    // compute relative tangential velocity
    // vector cVeliNorm = wallCntvar.Veli_
        // - ((wallCntvar.Veli_ & wallCntvar.contactNormal_)
        // *wallCntvar.contactNormal_);
    vector cVeliNorm = wallCntvar.Veli_*(wallCntvar.Veli_&wallCntvar.contactNormal_);

    vector Vt(wallCntvar.Veli_-(cVeliNorm - vector::zero));
    // compute tangential force
    if(contact_model_info::getUseMindlinRotationalModel())
    {
        
        scalar kT = 200*8*meanCntPar.aG_*(wallCntvar.contactArea_/(wallCntvar.Lc_+SMALL));
        vector deltaFt(kT*Vt*deltaT + 2*meanCntPar.reduceBeta_*sqrt(kT*reduceM_)*Vt);
        wallCntvar.FtPrev_ = - FtLastS - deltaFt;
    }

    if(contact_model_info::getUseChenRotationalModel())
    {
   
        vector Ftdi(meanCntPar.reduceBeta_*sqrt(meanCntPar.aG_*reduceM_*wallCntvar.Lc_)*Vt);
        Ftdi += meanCntPar.aG_*wallCntvar.Lc_*Vt*deltaT;
        wallCntvar.FtPrev_ = - FtLastS - Ftdi;
    }

    return wallCntvar.FtPrev_;
}
//---------------------------------------------------------------------------//
void wallSubContactInfo::syncData()
{
    reduce(outForce_.F, sumOp<vector>());
    reduce(outForce_.T, sumOp<vector>());
}
//---------------------------------------------------------------------------//
void wallSubContactInfo::syncContactResolve()
{
    reduce(contact_resolved_,orOp<bool>());
}
//---------------------------------------------------------------------------//
autoPtr<virtual_meshWallInfo>& wallSubContactInfo::getVMContactInfo
(
    label ID
)
{
    return vmWInfoList_[ID];
}
//---------------------------------------------------------------------------//
autoPtr<virtual_meshWallInfo>& wallSubContactInfo::getVMPlaneInfo
(
    label ID
)
{
    return vmPlaneInfoList_[ID];
}

//---------------------------------------------------------------------------//
