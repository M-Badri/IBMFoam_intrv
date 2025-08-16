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
#include "prtSubContactInfo.H"
#include "contact_model_info.H"

using namespace Foam;

//---------------------------------------------------------------------------//
prtSubContactInfo::prtSubContactInfo
(
    const Tuple2<label,label>& contactPair,
    const physicalProperties& physicalProperties
)
:
contactPair_(contactPair),
physicalProperties_(physicalProperties)
{}

prtSubContactInfo::~prtSubContactInfo()
{}
//---------------------------------------------------------------------------//
vector prtSubContactInfo::getVeli(ibContactVars& cVars, vector& lVec)
{
    return (-((lVec-cVars.Axis_*((lVec) & cVars.Axis_))
        ^ cVars.Axis_)*cVars.omega_+ cVars.Vel_);
}
//---------------------------------------------------------------------------//
void prtSubContactInfo::evalVariables(
    ibContactClass& cIb,
    ibContactClass& tIb,
    ibContactVars& cVars,
    ibContactVars& tVars
)
{
    cLVec_ = cIb.get_geom_model().getLVec(prtCntVars_.contactCenter_);
    tLVec_ = tIb.get_geom_model().getLVec(prtCntVars_.contactCenter_);
   

    cVeli_ = getVeli(cVars, cLVec_);
    tVeli_ = getVeli(tVars, tLVec_);

    Vn_ = -(cVeli_ - tVeli_) & prtCntVars_.contactNormal_;
    Lc_ = (contact_model_info::getLcCoeff())*mag(cLVec_)*mag(tLVec_)/(mag(cLVec_) + mag(tLVec_));

    physicalProperties_.curAdhN_ = min
    (
        physicalProperties_.maxAdhN_,
        max(physicalProperties_.curAdhN_, physicalProperties_.aY_*prtCntVars_.contactVolume_/(sqr(Lc_)*8
            *Foam::constant::mathematical::pi))
    );
}
//---------------------------------------------------------------------------//
vector prtSubContactInfo::getFNe()
{
    return (physicalProperties_.aY_*prtCntVars_.contactVolume_/(Lc_+SMALL))
        *prtCntVars_.contactNormal_;
}
//---------------------------------------------------------------------------//
vector prtSubContactInfo::getFA()
{
    return ((sqrt(8*Foam::constant::mathematical::pi*physicalProperties_.aY_
        *physicalProperties_.curAdhN_*prtCntVars_.contactVolume_))
        *prtCntVars_.contactNormal_);
}
//---------------------------------------------------------------------------//
vector prtSubContactInfo::getFNd()
{
    return (physicalProperties_.reduceBeta_*sqrt(physicalProperties_.aY_
        *physicalProperties_.reduceM_*prtCntVars_.contactArea_/(Lc_+SMALL))*
        Vn_)*prtCntVars_.contactNormal_;

}
//---------------------------------------------------------------------------//
vector prtSubContactInfo::getFt(scalar deltaT)
{
    vector FtLastP(FtPrev_ - (FtPrev_ & prtCntVars_.contactNormal_)
        *prtCntVars_.contactNormal_);


    vector FtLastS(mag(FtPrev_) * (FtLastP/(mag(FtLastP)+SMALL)));

    vector relVeli(cVeli_ - tVeli_);
    vector veliNomr((relVeli)*(relVeli & prtCntVars_.contactNormal_));
    vector Vt(relVeli-veliNomr);

    if(contact_model_info::getUseMindlinRotationalModel())
    {
        
        scalar kT = 200*8*physicalProperties_.aG_*(prtCntVars_.contactArea_/(Lc_+SMALL));
        vector deltaFt(kT*Vt*deltaT + 2*physicalProperties_.reduceBeta_*sqrt(kT*physicalProperties_.reduceM_)*Vt);
        FtPrev_ = - FtLastS - deltaFt;
    }

    if(contact_model_info::getUseChenRotationalModel())
    {
        
        vector Ftdi(- physicalProperties_.reduceBeta_*sqrt(physicalProperties_.aG_*physicalProperties_.reduceM_*Lc_)*Vt);
        Ftdi += physicalProperties_.aG_*Lc_*Vt*deltaT;
        FtPrev_ = - FtLastS- Ftdi;
    }

    
    return FtPrev_;
}
//---------------------------------------------------------------------------//
void prtSubContactInfo::setVMInfo(boundBox& bBox, scalar sub_volumeV)
{
    if (!vmInfo_)
    {
        vmInfo_ = std::make_shared<virtual_meshInfo>(bBox, sub_volumeV);
        return;
    }

    vmInfo_->sV = sub_volume(bBox);
    vmInfo_->sub_volumeV = sub_volumeV;
}
//---------------------------------------------------------------------------//
void prtSubContactInfo::setVMInfo(const virtual_meshInfo& vmInfo)
{
    if (!vmInfo_)
    {
        vmInfo_ = std::make_shared<virtual_meshInfo>(vmInfo);
        return;
    }

    vmInfo_->sV = vmInfo.sV;
    vmInfo_->sub_volumeV = vmInfo.sub_volumeV;
    vmInfo_->startingPoint = vmInfo.startingPoint;
}
//---------------------------------------------------------------------------//
std::shared_ptr<virtual_meshInfo>& prtSubContactInfo::getVMInfo()
{
    return vmInfo_;
}
//---------------------------------------------------------------------------//
void prtSubContactInfo::syncData()
{
    reduce(outForce_.first().F, sumOp<vector>());
    reduce(outForce_.first().T, sumOp<vector>());
    reduce(outForce_.second().F, sumOp<vector>());
    reduce(outForce_.second().T, sumOp<vector>());


    if (vmInfo_)
    {
        point reducePoint = vector::zero;

        if (contact_resolved_)
        {
            reducePoint = vmInfo_->getStartingPoint();
        }

        reduce(reducePoint, sumOp<vector>());
        vmInfo_->startingPoint.reset(new point(reducePoint));
    }
}
//---------------------------------------------------------------------------//


// ************************************************************************* //
