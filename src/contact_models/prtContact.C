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
#include "prtContact.H"

#include "virtual_meshLevel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace contactModel
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
//---------------------------------------------------------------------------//
void getContacts(
    const fvMesh&   mesh,
    prtContactInfo& prtcInfo
)
{
    if (
        prtcInfo.getcClass().get_geom_model().getcType() == sphere
        &&
        prtcInfo.gettClass().get_geom_model().getcType() == sphere
        )
    {
        prtcInfo.getContacts_Sphere();
    }
    else if (
        prtcInfo.getcClass().get_geom_model().getcType() == cluster
        ||
        prtcInfo.gettClass().get_geom_model().getcType() == cluster
        )
    {
        getContacts_Cluster(mesh, prtcInfo);
    }
    else
    {

        scalar ranCellVol = pow(virtual_meshLevel::getCharCellSize(),3);
        prtcInfo.getContacts_ArbShape(ranCellVol);        
    }

    prtcInfo.swapContactLists();
}
//---------------------------------------------------------------------------//
void getContacts_Cluster(
    const fvMesh&   mesh,
    prtContactInfo& prtcInfo
)
{
    prtcInfo.swapContactLists();
    std::vector<std::shared_ptr<geom_model>> cBodies;
    std::vector<std::shared_ptr<geom_model>> tBodies;

    bool cIsCluster = prtcInfo.getcClass().get_geom_model().isCluster();

    if(cIsCluster)
    {
        periodicBody& cCluster = dynamic_cast<periodicBody&>(prtcInfo.getcClass().get_geom_model());
        cBodies = cCluster.get_cluster_bodies();
    }
    else
    {
        cBodies.push_back(prtcInfo.getcClass().get_geom_modelPtr());
    }

    if(prtcInfo.gettClass().get_geom_model().isCluster())
    {
        periodicBody& tCluster = dynamic_cast<periodicBody&>(prtcInfo.gettClass().get_geom_model());

        if (cIsCluster)
        {
            tBodies.push_back(tCluster.get_cluster_bodies()[0]);
        }
        else
        {
            tBodies = tCluster.get_cluster_bodies();
        }
    }
    else
    {
        tBodies.push_back(prtcInfo.gettClass().get_geom_modelPtr());
    }

    for(std::shared_ptr<geom_model>& cgModel : cBodies)
    {
        for(std::shared_ptr<geom_model>& tgModel : tBodies)
        {
            ibContactClass cIbClassI(
                cgModel,
                prtcInfo.getcClass().getMatInfo().getMaterial()
            );

            ibContactClass tIbClassI(
                tgModel,
                prtcInfo.gettClass().getMatInfo().getMaterial()
            );

            prtContactInfo tmpPrtCntInfo
            (
                cIbClassI,
                prtcInfo.getcVars(),
                tIbClassI,
                prtcInfo.gettVars()
            );

            getContacts(mesh, tmpPrtCntInfo);

            if (tmpPrtCntInfo.getPrtSCList().size() > 0)
            {
                prtcInfo.getPrtSCList().insert(
                    prtcInfo.getPrtSCList().end(),
                    tmpPrtCntInfo.getPrtSCList().begin(),
                    tmpPrtCntInfo.getPrtSCList().end()
                );
            }
        }
    }
    prtcInfo.swapContactLists();
}
//---------------------------------------------------------------------------//
bool detectPrtPrtContact(
    const fvMesh&   mesh,
    ibContactClass& cClass,
    ibContactClass& tClass,
    prtSubContactInfo& subCInfo
)
{
    if
    (
        cClass.get_geom_model().getcType() == sphere
        &&
        tClass.get_geom_model().getcType() == sphere
    )
    {
        return detectPrtPrtContact_Sphere
        (
            mesh,
            cClass,
            tClass
        );
    }
    else if
    (
        cClass.get_geom_model().getcType() == cluster
        ||
        tClass.get_geom_model().getcType() == cluster
    )
    {
        return detectPrtPrtContact_Cluster
        (
            mesh,
            cClass,
            tClass,
            subCInfo
        );
    }
    else
    {
        return detectPrtPrtContact_ArbShape
        (
            mesh,
            cClass,
            tClass,
            subCInfo
        );
    }
}
//---------------------------------------------------------------------------//
bool detectPrtPrtContact_ArbShape(
    const fvMesh&   mesh,
    ibContactClass& cClass,
    ibContactClass& tClass,
    prtSubContactInfo& subCInfo
)
{
    std::shared_ptr<virtual_meshInfo>& vmInfo = subCInfo.getVMInfo();
    if (!vmInfo)
    {
        return false;
    }

    virtual_mesh virtMesh(
        *vmInfo,
        cClass.get_geom_model(),
        tClass.get_geom_model()
    );

    return virtMesh.detectFirstContactPoint();
}
//---------------------------------------------------------------------------//
bool detectPrtPrtContact_Sphere
(
    const fvMesh&   mesh,
    ibContactClass& cClass,
    ibContactClass& tClass
)
{
    if
    (
        mag(cClass.get_geom_model().getCoM()-tClass.get_geom_model().getCoM())
        <
        ((cClass.get_geom_model().getDC() / 2) + (tClass.get_geom_model().getDC() / 2))
    )
    {
        return true;
    }
    return false;
}
//---------------------------------------------------------------------------//
bool detectPrtPrtContact_Cluster
(
    const fvMesh&   mesh,
    ibContactClass& cClass,
    ibContactClass& tClass,
    prtSubContactInfo& subCInfo
)
{
    std::vector<std::shared_ptr<geom_model>> cBodies;
    std::vector<std::shared_ptr<geom_model>> tBodies;

    bool isCCluster = cClass.get_geom_model().isCluster();

    if(isCCluster)
    {
        periodicBody& cCluster = dynamic_cast<periodicBody&>(cClass.get_geom_model());
        cBodies = cCluster.get_cluster_bodies();
    }
    else
    {
        cBodies.push_back(cClass.get_geom_modelPtr());
    }

    if(tClass.get_geom_model().isCluster())
    {
        periodicBody& tCluster = dynamic_cast<periodicBody&>(tClass.get_geom_model());
        if (isCCluster)
        {
            tBodies.push_back(tCluster.get_cluster_bodies()[0]);
        }
        else
        {
            tBodies = tCluster.get_cluster_bodies();
        }
    }
    else
    {
        tBodies.push_back(tClass.get_geom_modelPtr());
    }

    for(std::shared_ptr<geom_model>& cgModel : cBodies)
    {
        for(std::shared_ptr<geom_model>& tgModel : tBodies)
        {
            ibContactClass cIbClassI(
                cgModel,
                cClass.getMatInfo().getMaterial()
            );

            ibContactClass tIbClassI(
                tgModel,
                tClass.getMatInfo().getMaterial()
            );

            prtSubContactInfo tmpSubCInfoI(
                subCInfo.getCPair(),
                subCInfo.getPhysicalProperties()
            );

            if (detectPrtPrtContact(
                mesh,
                cIbClassI,
                tIbClassI,
                subCInfo
            ))
            {
                return true;
            }
        }
    }
    return false;
}
//---------------------------------------------------------------------------//
void getPrtContactVars_ArbShape(
    const fvMesh&   mesh,
    ibContactClass& cClass,
    ibContactClass& tClass,
    prtSubContactInfo& subCInfo
)
{
    std::shared_ptr<virtual_meshInfo>& vmInfo = subCInfo.getVMInfo();
    if (!vmInfo)
    {
        return;
    }

    virtual_mesh virtMesh(
        *vmInfo,
        cClass.get_geom_model(),
        tClass.get_geom_model()
    );

    scalar intersectedVolume(0);
    vector contactCenter(vector::zero);
    vector normalVector = vector::zero;
    scalar contactArea(0);

    intersectedVolume = virtMesh.evaluateContact();

    if(virtMesh.getEdgeSVPoints().size() <= 4)
    {
        intersectedVolume = 0;
    }

    if (intersectedVolume > 0)
    {
        Tuple2<scalar,vector> surfaceAndNormal =
            virtMesh.get3DcontactNormalAndSurface
            (
                cClass.get_geom_model().getcType() == nonConvex
                ||
                tClass.get_geom_model().getcType() == nonConvex
            );
        contactArea = surfaceAndNormal.first();
        normalVector = surfaceAndNormal.second();
        contactCenter = virtMesh.getContactCenter();
    }

    if (intersectedVolume > 0 && contactArea > 0)
    {
        subCInfo.getprtCntVars().contactCenter_ = contactCenter;
        subCInfo.getprtCntVars().contactVolume_ = intersectedVolume;
        subCInfo.getprtCntVars().contactNormal_ = normalVector;
        subCInfo.getprtCntVars().contactArea_ = contactArea;
    }
    else
    {
        subCInfo.getprtCntVars().contactCenter_ = vector::zero;
        subCInfo.getprtCntVars().contactVolume_ = 0;
        subCInfo.getprtCntVars().contactNormal_ = vector::zero;
        subCInfo.getprtCntVars().contactArea_ = 0;
    }
}
//---------------------------------------------------------------------------//
void getPrtContactVars_Sphere
(
    const fvMesh&   mesh,
    ibContactClass& cClass,
    ibContactClass& tClass,
    prtSubContactInfo& subCInfo
)
{
    scalar cRadius = cClass.get_geom_model().getDC() / 2;
    scalar tRadius = tClass.get_geom_model().getDC() / 2;
    vector tCenter = tClass.get_geom_model().getCoM();

    vector centerDir = cClass.get_geom_model().getCoM()
                        - tCenter;

    scalar d = mag(centerDir);

    if(mag(centerDir) < SMALL || d > (cRadius + tRadius))
    {
        subCInfo.getprtCntVars().contactCenter_ = vector::zero;
        subCInfo.getprtCntVars().contactVolume_ = 0;
        subCInfo.getprtCntVars().contactNormal_ = vector::zero;
        subCInfo.getprtCntVars().contactArea_ = 0;
        return;
    }

    scalar xLength = (sqr(d) - sqr(cRadius) + sqr(tRadius))/(2*d);

    if(sqr(xLength) > sqr(tRadius))
    {
        subCInfo.getprtCntVars().contactCenter_
            = tCenter + (centerDir/d)*xLength;
        subCInfo.getprtCntVars().contactVolume_
            = (4/3)*Foam::constant::mathematical::pi*pow(tRadius,3);
        subCInfo.getprtCntVars().contactNormal_ = centerDir/d;
        subCInfo.getprtCntVars().contactArea_
            = Foam::constant::mathematical::pi*sqr(tRadius);
        return;
    }

    if(case3D)
    {

        scalar cSphCapV = (Foam::constant::mathematical::pi
                            *sqr(cRadius - d + xLength)
                            *(3*cRadius - (cRadius - d + xLength))) / 3;
        scalar tSphCapV = (Foam::constant::mathematical::pi
                            *sqr(tRadius - xLength)
                            *(3*tRadius - (tRadius - xLength))) / 3;

        subCInfo.getprtCntVars().contactCenter_
            = tCenter + (centerDir/d)*xLength;
        subCInfo.getprtCntVars().contactVolume_
            = cSphCapV + tSphCapV;
        subCInfo.getprtCntVars().contactNormal_
            = centerDir/d;
        subCInfo.getprtCntVars().contactArea_
            = Foam::constant::mathematical::pi
                *sqr(sqrt(sqr(tRadius) - sqr(xLength)));
    }
    else
    {
        boundBox meshBounds = mesh.bounds();
        scalar emptyLength = meshBounds.max()[emptyDim]
                            - meshBounds.min()[emptyDim];

        scalar cCirSeg = sqr(cRadius)*acos((d - xLength)/cRadius)
                         - (d - xLength)*sqrt(sqr(cRadius)
                         - sqr(d - xLength));
        scalar tCirSeg = sqr(tRadius)*acos(xLength/tRadius)
                         - xLength*sqrt(sqr(tRadius) - sqr(xLength));

        subCInfo.getprtCntVars().contactCenter_
            = tCenter + (centerDir/d)*xLength;
        subCInfo.getprtCntVars().contactVolume_
            = (cCirSeg + tCirSeg)*emptyLength;
        subCInfo.getprtCntVars().contactNormal_
            = centerDir/d;
        subCInfo.getprtCntVars().contactArea_
            = 2*sqrt(sqr(tRadius) - sqr(xLength))*emptyLength;
    }
}
//---------------------------------------------------------------------------//
void getPrtContactVars_Cluster
(
    const fvMesh&   mesh,
    ibContactClass& cClass,
    ibContactClass& tClass,
    prtSubContactInfo& subCInfo
)
{
    subCInfo.getprtCntVars().contactCenter_ = vector::zero;
    subCInfo.getprtCntVars().contactVolume_ = 0;
    subCInfo.getprtCntVars().contactNormal_ = vector::zero;
    subCInfo.getprtCntVars().contactArea_ = 0;

    std::vector<std::shared_ptr<geom_model>> cBodies;
    std::vector<std::shared_ptr<geom_model>> tBodies;

    bool isCCluster = cClass.get_geom_model().isCluster();

    if(isCCluster)
    {
        periodicBody& cCluster = dynamic_cast<periodicBody&>(cClass.get_geom_model());
        cBodies = cCluster.get_cluster_bodies();
    }
    else
    {
        cBodies.push_back(cClass.get_geom_modelPtr());
    }

    if(tClass.get_geom_model().isCluster())
    {
        periodicBody& tCluster = dynamic_cast<periodicBody&>(tClass.get_geom_model());
        if (isCCluster)
        {
            tBodies.push_back(tCluster.get_cluster_bodies()[0]);
        }
        else
        {
            tBodies = tCluster.get_cluster_bodies();
        }
    }
    else
    {
        tBodies.push_back(tClass.get_geom_modelPtr());
    }

    for(std::shared_ptr<geom_model>& cgModel : cBodies)
    {
        for(std::shared_ptr<geom_model>& tgModel : tBodies)
        {
            ibContactClass cIbClassI(
                cgModel,
                cClass.getMatInfo().getMaterial()
            );

            ibContactClass tIbClassI(
                tgModel,
                tClass.getMatInfo().getMaterial()
            );

            prtSubContactInfo tmpSubCInfoI(
                subCInfo.getCPair(),
                subCInfo.getPhysicalProperties()
            );

            if (subCInfo.getVMInfo())
            {
                tmpSubCInfoI.setVMInfo(*(subCInfo.getVMInfo()));
            }

            getPrtContactVars(
                mesh,
                cIbClassI,
                tIbClassI,
                tmpSubCInfoI
            );

            if (tmpSubCInfoI.getprtCntVars().contactVolume_
                > subCInfo.getprtCntVars().contactVolume_)
            {
                subCInfo.getprtCntVars() = tmpSubCInfoI.getprtCntVars();
            }
        }
    }
}
//---------------------------------------------------------------------------//
void getPrtContactVars
(
    const fvMesh&   mesh,
    ibContactClass& cClass,
    ibContactClass& tClass,
    prtSubContactInfo& subCInfo
)
{
    if
    (
        cClass.get_geom_model().getcType() == sphere
        &&
        tClass.get_geom_model().getcType() == sphere
    )
    {
        getPrtContactVars_Sphere(
            mesh,
            cClass,
            tClass,
            subCInfo
        );
    }
    else if
    (
        cClass.get_geom_model().getcType() == cluster
        ||
        tClass.get_geom_model().getcType() == cluster
    )
    {
        getPrtContactVars_Cluster(
            mesh,
            cClass,
            tClass,
            subCInfo
        );
    }
    else
    {
        getPrtContactVars_ArbShape(
            mesh,
            cClass,
            tClass,
            subCInfo
        );
    }
}
//---------------------------------------------------------------------------//
bool solvePrtContact(
    const fvMesh&   mesh,
    prtContactInfo& cInfo,
    prtSubContactInfo& subCInfo,
    scalar deltaT
)
{
    getPrtContactVars(
        mesh,
        cInfo.getcClass(),
        cInfo.gettClass(),
        subCInfo
    );

    if (!(subCInfo.getprtCntVars().contactVolume_ > 0))
    {
        return false;
    }

    InfoH << parallelDEM_Info << "-- Detected Particle-particle contact: -- body "
            << subCInfo.getCPair().first() << " & -- body "
            << subCInfo.getCPair().second() << endl;
    InfoH << parallelDEM_Info << "-- Particle-particle " <<subCInfo.getCPair().first() <<"-"<<subCInfo.getCPair().second() << " contact cBody pos: "
            << cInfo.getcClass().get_geom_model().getCoM() << " & tBody pos: "
            << cInfo.gettClass().get_geom_model().getCoM() << endl;
    InfoH << parallelDEM_Info << "-- body "<< subCInfo.getCPair().first() <<"  linear velocity:"
        << cInfo.getcVars().Vel_ << " magnitude: " << mag(cInfo.getcVars().Vel_) <<endl;
    InfoH << parallelDEM_Info << "-- body "<< subCInfo.getCPair().first() <<"  angular velocity:"
        << cInfo.getcVars().omega_ << " magnitude: " << mag(cInfo.getcVars().omega_) <<endl;

    InfoH << parallelDEM_Info << "-- body "<< subCInfo.getCPair().second() <<"  linear velocity:"
        << cInfo.gettVars().Vel_ << " magnitude: " << mag(cInfo.gettVars().Vel_) <<endl;
    InfoH << parallelDEM_Info << "-- body "<< subCInfo.getCPair().second() <<"  angular velocity:"
        << cInfo.gettVars().omega_ << " magnitude: " << mag(cInfo.gettVars().omega_) <<endl;
    InfoH << parallelDEM_Info << "-- Particle-particle " <<subCInfo.getCPair().first() <<"-"<<subCInfo.getCPair().second() << " contact center "
            << subCInfo.getprtCntVars().contactCenter_ << endl;
    InfoH << parallelDEM_Info << "-- Particle-particle " <<subCInfo.getCPair().first() <<"-"<<subCInfo.getCPair().second() << " contact normal "
            << subCInfo.getprtCntVars().contactNormal_ << endl;
    InfoH << parallelDEM_Info << "-- Particle-particle " <<subCInfo.getCPair().first() <<"-"<<subCInfo.getCPair().second() << " contact volume "
            << subCInfo.getprtCntVars().contactVolume_ << endl;
    InfoH << parallelDEM_Info << "-- Particle-particle " <<subCInfo.getCPair().first() <<"-"<<subCInfo.getCPair().second() << " contact area "
            << subCInfo.getprtCntVars().contactArea_ << endl;

    subCInfo.evalVariables(
        cInfo.getcClass(),
        cInfo.gettClass(),
        cInfo.getcVars(),
        cInfo.gettVars()
    );

    // compute the normal force
    vector F = subCInfo.getFNe();
    InfoH << parallelDEM_Info << "-- Particle-particle " <<subCInfo.getCPair().first() <<"-"<<subCInfo.getCPair().second() << " contact FNe " << F << endl;

    vector FNd = subCInfo.getFNd();
    InfoH << parallelDEM_Info << "-- Particle-particle " <<subCInfo.getCPair().first() <<"-"<<subCInfo.getCPair().second() << " contact FNd " << FNd << endl;

    // clamp FNd if opposite direction to FNe
    if ((F & FNd) < 0 && mag(FNd) > mag(F))
    {
        FNd *= mag(F) / mag(FNd);
    }

    F += FNd;
    InfoH << parallelDEM_Info << "-- Particle-particle " <<subCInfo.getCPair().first() <<"-"<<subCInfo.getCPair().second() << " contact FN " << F << endl;

    vector Ft = subCInfo.getFt(deltaT);
    InfoH << parallelDEM_Info << "-- Particle-particle " <<subCInfo.getCPair().first() <<"-"<<subCInfo.getCPair().second() << " contact Ft " << Ft << endl;

    if (mag(Ft) > cInfo.getMu() * mag(F))
    {
        Ft *= cInfo.getMu() * mag(F) / mag(Ft);
    }
    InfoH << parallelDEM_Info << "-- Particle-particle " <<subCInfo.getCPair().first() <<"-"<<subCInfo.getCPair().second() << " contact Ft clamped " << Ft << endl;
    F += Ft;

    vector FA = subCInfo.getFA();
    InfoH << parallelDEM_Info << "-- Particle-particle " <<subCInfo.getCPair().first() <<"-"<<subCInfo.getCPair().second() << " contact FA " << FA << endl;
    F -= FA;


    subCInfo.get_out_force().first().F = F;
    subCInfo.get_out_force().first().T = subCInfo.getcLVec() ^  F;
    subCInfo.get_out_force().second().F = -F;
    subCInfo.get_out_force().second().T = subCInfo.gettLVec() ^ -F;

    InfoH << parallelDEM_Info << "-- Particle-particle " <<subCInfo.getCPair().first() <<"-"<<subCInfo.getCPair().second() << " contact F " << F << endl;

    InfoH << parallelDEM_Info << "-- Particle-particle " <<subCInfo.getCPair().first() <<"-"<<subCInfo.getCPair().second() << " contact T " << subCInfo.get_out_force().first().T << endl;

    InfoH << parallelDEM_Info << "-- Resolved Particle-particle contact: -- body "
            << subCInfo.getCPair().first() << " & -- body "
            << subCInfo.getCPair().second() << endl;

    return true;
}
//---------------------------------------------------------------------------//
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace contactModel

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
