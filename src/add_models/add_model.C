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
#include "add_model.H"

using namespace Foam;

//---------------------------------------------------------------------------//
add_model::~add_model()
{
}
//---------------------------------------------------------------------------//
bool add_model::isBodyInContact(PtrList<immersed_body>& imm_bodies)
{
    std::shared_ptr<geom_model> gModel = geom_model_->getCopy();
    ibContactClass cIbClass(
        gModel,
        "None"
    );

    vector velAxi (vector::zero);
    scalar helpScalar(0.0);

    ibContactVars cIbVars(
        imm_bodies.size(),
        velAxi,
        helpScalar,
        velAxi,
        helpScalar,
        helpScalar,
        geom_model_->getRhoS()
    );

    wallContactInfo cIBWallCntI(
            cIbClass,
            cIbVars
    );

    bool inContact = contactModel::detectWallContact(
        mesh_,
        cIbClass,
        cIBWallCntI
    );

    if(!inContact)
    {
        forAll(imm_bodies, ibI)
        {
            bool bBoxContact = true;
            boundBox ibIbBox = imm_bodies[ibI].get_geom_model().getBounds();
            boundBox cbBox = geom_model_->getBounds();

            forAll(geometricD,dir)
            {
                if(geometricD[dir] == 1)
                {
                    if(!(ibIbBox.max()[dir] >= cbBox.min()[dir]
                        && ibIbBox.min()[dir] <= cbBox.max()[dir]))
                    {
                        bBoxContact = false;
                        break;
                    }
                }
            }

            if(!bBoxContact)
            {
                continue;
            }

            prtContactInfo prtCInfo (
                cIbClass,
                cIbVars,
                imm_bodies[ibI].get_ib_contact_class(),
                imm_bodies[ibI].getContactVars()
            );

            contactModel::getContacts(
                mesh_,
                prtCInfo
            );

            DynamicList<prtSubContactInfo*> sCList;
            prtCInfo.registerContactList(sCList);
            forAll(sCList,sC)
            {
                prtSubContactInfo* prtSCInfo = sCList[sC];

                if(contactModel::detectPrtPrtContact(
                    mesh_,
                    cIbClass,
                    imm_bodies[ibI].get_ib_contact_class(),
                    *prtSCInfo
                ))
                {
                    inContact = true;
                    break;
                }

            }

        }
    }

    return inContact;
}
