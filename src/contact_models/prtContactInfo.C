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
#include "prtContactInfo.H"

#include "inter_adhesion.H"
#include "virtual_meshLevel.H"
#include "contact_model_info.H"

using namespace Foam;

//---------------------------------------------------------------------------//
prtContactInfo::prtContactInfo
(
    ibContactClass& cClass,
    ibContactVars& cVars,
    ibContactClass& tClass,
    ibContactVars& tVars
)
:
cIbContactClass_(cClass),
cContactVars_(cVars),
tIbContactClass_(tClass),
tContactVars_(tVars)
{
    contactPair_.first() = cVars.body_id_;
    contactPair_.second() = tVars.body_id_;
    const material_info& cMatInfo(cClass.getMatInfo());
    const material_info& tMatInfo(tClass.getMatInfo());

    scalar adhPot = 0;
    string adhPotKey;
    if(cMatInfo.getMaterial() < tMatInfo.getMaterial())
    {
        adhPotKey = cMatInfo.getMaterial();
        adhPotKey += "-";
        adhPotKey += tMatInfo.getMaterial();
    }
    else
    {
        adhPotKey = tMatInfo.getMaterial();
        adhPotKey += "-";
        adhPotKey += cMatInfo.getMaterial();
    }

    if(inter_adhesion::getInterAdhesion().found(adhPotKey))
    {
        adhPot = inter_adhesion::getInterAdhesion()[adhPotKey];
    }

    // compute mean model parameters
    physicalProperties_.aY_ = 1/((1 - sqr(cMatInfo.getNu()))/cMatInfo.getY()
        + (1 - sqr(tMatInfo.getNu()))/tMatInfo.getY());
    physicalProperties_.aG_ = 1/(2*(2 - cMatInfo.getNu())*(1 + cMatInfo.getNu())/cMatInfo.getY()
        + 2*(2 - tMatInfo.getNu())*(1 + tMatInfo.getNu())/tMatInfo.getY());
    physicalProperties_.aMu_ = (cMatInfo.getMu()+tMatInfo.getMu())/2;
    physicalProperties_.maxAdhN_ = cMatInfo.getAdhN() + tMatInfo.getAdhN() - 2*adhPot;
    physicalProperties_.curAdhN_ = 0;
    physicalProperties_.reduceM_ =
    (
        cIbContactClass_.get_geom_model().getM0()
        *tIbContactClass_.get_geom_model().getM0()
        /(cIbContactClass_.get_geom_model().getM0()
        +tIbContactClass_.get_geom_model().getM0())
    );
    
    physicalProperties_.reduceBeta_ =
    (
       (-1)*sqrt(5.0)*log((0.5*(cMatInfo.getEps()+tMatInfo.getEps())))/
       (sqrt(sqr(log(0.5*(cMatInfo.getEps()+tMatInfo.getEps())))+
       sqr(Foam::constant::mathematical::pi)))  
    );
}

prtContactInfo::~prtContactInfo()
{
}
//---------------------------------------------------------------------------//
std::shared_ptr<prtSubContactInfo> prtContactInfo::matchSubContact
(
    boundBox& bbox,
    physicalProperties& physicalProperties,
    Tuple2<label,label>& contactPair
)
{
    for(auto sC : contactList_)
    {
        if (!sC->getVMInfo())
        {
            continue;
        }

        if (bbox.contains(sC->getVMInfo()->getStartingPoint()))
        {
            return sC;
        }
    }

    return std::make_shared<prtSubContactInfo>
        (contactPair, physicalProperties);
}
//---------------------------------------------------------------------------//
void prtContactInfo::limitBBox(boundBox& bbox)
{
    const boundBox cBBox(cIbContactClass_.get_geom_model().getBounds());
    const boundBox tBBox(tIbContactClass_.get_geom_model().getBounds());
    for (label coord = 0; coord < 3; coord++)
    {
        if (bbox.min()[coord] < cBBox.min()[coord])
        {
            bbox.min()[coord] = cBBox.min()[coord];
        }

        if (bbox.max()[coord] > cBBox.max()[coord])
        {
            bbox.max()[coord] = cBBox.max()[coord];
        }

        if (bbox.min()[coord] < tBBox.min()[coord])
        {
            bbox.min()[coord] = tBBox.min()[coord];
        }

        if (bbox.max()[coord] > tBBox.max()[coord])
        {
            bbox.max()[coord] = tBBox.max()[coord];
        }
    }
}
//---------------------------------------------------------------------------//
void prtContactInfo::getContacts_Sphere()
{
    if
    (
        mag(cIbContactClass_.get_geom_model().getCoM()-tIbContactClass_.get_geom_model().getCoM())
        >=
        ((cIbContactClass_.get_geom_model().getDC() / 2) + (tIbContactClass_.get_geom_model().getDC() / 2))
    )
    {
        return;
    }

    newContactList_.emplace_back(std::make_shared<prtSubContactInfo>
        (contactPair_, physicalProperties_)
    );
}
//---------------------------------------------------------------------------//
void prtContactInfo::getContacts_ArbShape
(
    scalar cellV
)
{
    boundBox subCbBox(
        cIbContactClass_.get_geom_model().getBounds().min(),
        cIbContactClass_.get_geom_model().getBounds().max()
    );

    if (!subCbBox.overlaps(tIbContactClass_.get_geom_model().getBounds()))
    {
        return;
    }

    limitBBox(subCbBox);

    scalar charCellSize = pow(cellV,0.333333);
    scalar sub_volumeLength = charCellSize/virtual_meshLevel::getLevelOfDivision();
    scalar sub_volumeV = pow(sub_volumeLength,3);

    newContactList_.emplace_back(matchSubContact(subCbBox, physicalProperties_, contactPair_));
    newContactList_.back()->setVMInfo(subCbBox, sub_volumeV);
    return;
}
//---------------------------------------------------------------------------//
void prtContactInfo::swapContactLists()
{
    contactList_.swap(newContactList_);

}
//---------------------------------------------------------------------------//
bool prtContactInfo::contact_resolved()
{
    for(auto sC : contactList_)
    {
        if (sC->contact_resolved())
        {
            return true;
        }
    }
    return false;
}
//---------------------------------------------------------------------------//
void prtContactInfo::syncContactList()
{
    std::vector<std::shared_ptr<prtSubContactInfo>> syncedContactList;
    
    for (label procI = 0; procI < Pstream::nProcs(); procI++)
    {
        label numOfCnts = 0;
        if (procI == Pstream::myProcNo())
        {
            numOfCnts = contactList_.size();
        }
        reduce(numOfCnts, sumOp<label>());

        for (label i = 0; i < numOfCnts; i++)
        {
            bool vmInfoValid = false;
            if (procI == Pstream::myProcNo())
            {
                vmInfoValid = contactList_[i]->getVMInfo() ? true : false;
            }
            reduce(vmInfoValid, orOp<bool>());

            if (vmInfoValid)
            {
                syncedContactList.emplace_back(std::make_shared<prtSubContactInfo>
                    (contactPair_, physicalProperties_)
                );

                virtual_meshInfo vmInfoToSync;
                if (procI == Pstream::myProcNo())
                {
                    vmInfoToSync = virtual_meshInfo(*(contactList_[i]->getVMInfo()));
                }
                reduce(vmInfoToSync.sV.min(), sumOp<vector>());
                reduce(vmInfoToSync.sV.max(), sumOp<vector>());
                reduce(vmInfoToSync.sub_volumeV, sumOp<scalar>());

                point startPointToReduce = vmInfoToSync.getStartingPoint();
                if (procI != Pstream::myProcNo())
                {
                    startPointToReduce = vector::zero;
                }
                reduce(startPointToReduce, sumOp<vector>());

                syncedContactList.back()->setVMInfo(vmInfoToSync);
            }
            else if (procI == 0)
            {
                syncedContactList.emplace_back(std::make_shared<prtSubContactInfo>
                    (contactPair_, physicalProperties_)
                );
            }
        }
    }

    contactList_.swap(syncedContactList);
}
//---------------------------------------------------------------------------//
void prtContactInfo::registerContactList(DynamicList<prtSubContactInfo*>& contactList)
{
    for(auto sC : contactList_)
    {
        contactList.append(sC.get());
    }
}
//---------------------------------------------------------------------------//
void prtContactInfo::clearData()
{
    newContactList_.clear();
    for(auto sC : contactList_)
    {
        sC->clearOutForces();
        sC->setResolvedContact(false);
    }
}
//---------------------------------------------------------------------------//
void prtContactInfo::syncData()
{
    for(auto sC : contactList_)
    {
        sC->syncData();
    }
}
//---------------------------------------------------------------------------//

// ************************************************************************* //
