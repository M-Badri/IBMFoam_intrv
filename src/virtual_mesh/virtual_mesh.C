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
#include "virtual_mesh.H"

#include "sub_volume.H"
#include "virtual_meshLevel.H"

using namespace Foam;

//---------------------------------------------------------------------------//
virtual_mesh::virtual_mesh
(
    virtual_meshInfo& vMeshInfo,
    geom_model& cGeomModel,
    geom_model& tGeomModel
)
:
cGeomModel_(cGeomModel),
tGeomModel_(tGeomModel),
vMeshInfo_(vMeshInfo)
{}

virtual_mesh::~virtual_mesh()
{}
//---------------------------------------------------------------------------//
bool virtual_mesh::detectFirstContactPoint()
{
    bool startPointFound = false;

    return detectFirstVolumeInContact(
        vMeshInfo_.getSubVolume(),
        startPointFound
    );
}
//---------------------------------------------------------------------------//
bool virtual_mesh::detectFirstVolumeInContact(sub_volume& sV, bool& startPointFound)
{
    if (sV.volume() < vMeshInfo_.sub_volumeV)
    {
        if (sV.cVolumeInfo().volumeType_ != volumeType::inside)
        {
            sV.cVolumeInfo().volumeType_ = cGeomModel_.getVolumeType(sV, true);
        }

        if (sV.tVolumeInfo().volumeType_ != volumeType::inside)
        {
            sV.tVolumeInfo().volumeType_ = tGeomModel_.getVolumeType(sV, false);
        }

        if (sV.cVolumeInfo().volumeType_ == volumeType::outside
            || sV.tVolumeInfo().volumeType_ == volumeType::outside)
        {
            return false;
        }

        boundBox cBBox = boundBox(sV.min(), sV.max());
        boundBox tBBox = boundBox(sV.min(), sV.max());

        if (sV.cVolumeInfo().volumeType_ == volumeType::mixed)
        {
            if (!cGeomModel_.limitFinalSubVolume(
                sV,
                true,
                cBBox
            ))
            {
                return false;
            }
        }

        if (sV.tVolumeInfo().volumeType_ == volumeType::mixed)
        {
            if (!tGeomModel_.limitFinalSubVolume(
                sV,
                false,
                tBBox
            ))
            {
                return false;
            }
        }

        if (cBBox.overlaps(tBBox))
        {
            return true;
        }
        return false;
    }

    if (sV.cVolumeInfo().volumeType_ != volumeType::inside)
    {
        sV.cVolumeInfo().volumeType_ = cGeomModel_.getVolumeType(sV, true);
    }

    if (sV.tVolumeInfo().volumeType_ != volumeType::inside)
    {
        sV.tVolumeInfo().volumeType_ = tGeomModel_.getVolumeType(sV, false);
    }

    if (sV.cVolumeInfo().volumeType_ == volumeType::outside || sV.tVolumeInfo().volumeType_ == volumeType::outside)
    {
        return false;
    }
    else if (sV.cVolumeInfo().volumeType_ == volumeType::inside && sV.tVolumeInfo().volumeType_ == volumeType::inside)
    {
        return true;
    }

    List<sub_volume>& sVs = sV.childSubVolumes();

    if (!startPointFound)
    {
        for (auto iter = sVs.begin(); iter != sVs.end(); ++iter)
        {
            if (iter->contains(vMeshInfo_.getStartingPoint()))
            {
                // Move this one to the front
                std::rotate(sVs.begin(), iter, sVs.end());
                break;
            }
        }
    }

    forAll(sVs,i)
    {
        if (detectFirstVolumeInContact(sVs[i], startPointFound))
        {
            return true;
        }
        startPointFound = true;
    }
    return false;
}
//---------------------------------------------------------------------------//
scalar virtual_mesh::evaluateContact()
{
    scalar contactVolume = 0;
    edgeSubVolumesPoints_.clear();
    contactCenter_ = vector::zero;

    inspectSubVolume(vMeshInfo_.getSubVolume(), contactVolume, contactCenter_, edgeSubVolumesPoints_);

    if (contactVolume < VSMALL)
    {
        return 0;
    }

    contactCenter_ /= contactVolume;
    return contactVolume;
}
//---------------------------------------------------------------------------//
void virtual_mesh::inspectSubVolume(
    sub_volume& sV,
    scalar& contactVolume,
    vector& contactCenter,
    DynamicPointList& edgePoints
)
{
    if (sV.volume() < vMeshInfo_.sub_volumeV)
    {
        if (sV.cVolumeInfo().volumeType_ != volumeType::inside)
        {
            sV.cVolumeInfo().volumeType_ = cGeomModel_.getVolumeType(sV, true);
        }

        if (sV.tVolumeInfo().volumeType_ != volumeType::inside)
        {
            sV.tVolumeInfo().volumeType_ = tGeomModel_.getVolumeType(sV, false);
        }

        boundBox cBBox = boundBox(sV.min(), sV.max());
        boundBox tBBox = boundBox(sV.min(), sV.max());

        bool cVolumeTypeMixed = false;
        if (sV.cVolumeInfo().volumeType_ == volumeType::mixed)
        {
            cVolumeTypeMixed = true;

            sV.cVolumeInfo().volumeType_ = cGeomModel_.limitFinalSubVolume(
                sV,
                true,
                cBBox
            ) ? volumeType::inside : volumeType::outside;
        }

        if (sV.tVolumeInfo().volumeType_ == volumeType::mixed)
        {
            if (cVolumeTypeMixed)
            {
                edgePoints.append(sV.midpoint());
                sV.setAsEdge();
            }

            sV.tVolumeInfo().volumeType_ = tGeomModel_.limitFinalSubVolume(
                sV,
                false,
                tBBox
            ) ? volumeType::inside : volumeType::outside;
        }

        if (sV.cVolumeInfo().volumeType_ == volumeType::inside
            && sV.tVolumeInfo().volumeType_ == volumeType::inside
            && cBBox.overlaps(tBBox))
        {
            // boundBox as a cross section of cBBox and tBBox
            boundBox crossBBox(
                vector(
                    max(cBBox.min().x(), tBBox.min().x()),
                    max(cBBox.min().y(), tBBox.min().y()),
                    max(cBBox.min().z(), tBBox.min().z())
                ),
                vector(
                    min(cBBox.max().x(), tBBox.max().x()),
                    min(cBBox.max().y(), tBBox.max().y()),
                    min(cBBox.max().z(), tBBox.max().z())
                )
            );

            if (crossBBox.volume() > 0)
            {
                contactVolume += crossBBox.volume();
                contactCenter += (crossBBox.midpoint()*crossBBox.volume());
            }
        }
        return;
    }

    if (sV.cVolumeInfo().volumeType_ != volumeType::inside)
    {
        sV.cVolumeInfo().volumeType_ = cGeomModel_.getVolumeType(sV, true);
    }

    if (sV.tVolumeInfo().volumeType_ != volumeType::inside)
    {
        sV.tVolumeInfo().volumeType_ = tGeomModel_.getVolumeType(sV, false);
    }

    if (sV.cVolumeInfo().volumeType_ == volumeType::outside || sV.tVolumeInfo().volumeType_ == volumeType::outside)
    {
        return;
    }
    else if (sV.cVolumeInfo().volumeType_ == volumeType::inside && sV.tVolumeInfo().volumeType_ == volumeType::inside)
    {
        contactVolume += sV.volume();
        contactCenter += (sV.midpoint()*sV.volume());
        return;
    }

    List<sub_volume>& sVs = sV.childSubVolumes();
    DynamicPointList newEdgePoints;
    forAll(sVs,i)
    {
        inspectSubVolume(sVs[i], contactVolume, contactCenter, newEdgePoints);
    }

    if (newEdgePoints.size() > 0)
    {
        edgePoints.append(newEdgePoints);
        return;
    }

    if (sV.cVolumeInfo().volumeType_ == volumeType::mixed && sV.tVolumeInfo().volumeType_ == volumeType::mixed)
    {
        edgePoints.append(sV.midpoint());
        sV.setAsEdge();
    }
}
//---------------------------------------------------------------------------//
Tuple2<scalar,vector> virtual_mesh::get3DcontactNormalAndSurface(bool nonConvex)
{
    if (nonConvex)
    {
        std::vector<sub_contact> sCS = findsub_contacts(vMeshInfo_.sV);
        Tuple2<scalar,vector> contactNormalAndSurface
            = Tuple2<scalar,vector>(0.0,vector::zero);

        scalar totalVolume = 0;

        for (sub_contact& sC : sCS)
        {
            DynamicPointList edgeSubContactPoints = sC.getEdgePoints();
            if(edgeSubContactPoints.size() < 1)
            {
                continue;
            }
            Tuple2<scalar,vector> cNormalAndSurf
                = get3DcontactNormalAndSurface(edgeSubContactPoints);

            if (cNormalAndSurf.first() > 0)
            {
                contactNormalAndSurface.first() += cNormalAndSurf.first();
                contactNormalAndSurface.second()
                    += cNormalAndSurf.second() * sC.getVolume();
                totalVolume += sC.getVolume();
            }
        }

        if (contactNormalAndSurface.first() > 0)
        {
            contactNormalAndSurface.second()
                /= totalVolume;
        }

        return contactNormalAndSurface;
    }
    else
    {
        return get3DcontactNormalAndSurface(edgeSubVolumesPoints_);
    }
}
//---------------------------------------------------------------------------//
Tuple2<scalar,vector> virtual_mesh::get3DcontactNormalAndSurface(DynamicPointList edgeSubVolumesPoints)
{
// This function is taken from prtContact and just adjustated for higher accuracy
    scalar area(0.0);
    vector normalVec(vector::zero);
    scalar tDC(tGeomModel_.getDC());
    vector normalVector = vector::zero;
    point closestPoint = vector::zero;

    point sub_contactCenter = vector::zero;
    forAll (edgeSubVolumesPoints,cell)
    {
        sub_contactCenter += edgeSubVolumesPoints[cell];
    }
    sub_contactCenter /= edgeSubVolumesPoints.size();

    tGeomModel_.getClosestPointAndNormal
    (
        sub_contactCenter,
        vector::one * tDC,
        closestPoint,
        normalVector
    );

    if (edgeSubVolumesPoints.size() >= 3)
    {
        bool normOk(false);
        vector center(vector::zero);

        forAll (edgeSubVolumesPoints,cell)
        {
            center += edgeSubVolumesPoints[cell];
        }
        center /= edgeSubVolumesPoints.size();

        scalar xx(0);
        scalar xy(0);
        scalar xz(0);
        scalar yy(0);
        scalar yz(0);
        scalar zz(0);

        forAll (edgeSubVolumesPoints,cell)
        {
            vector subPoint(edgeSubVolumesPoints[cell] - center);
            if(subPoint != vector::zero)
                subPoint = subPoint/mag(subPoint);
            xx += subPoint[0] * subPoint[0];
            xy += subPoint[0] * subPoint[1];
            xz += subPoint[0] * subPoint[2];
            yy += subPoint[1] * subPoint[1];
            yz += subPoint[1] * subPoint[2];
            zz += subPoint[2] * subPoint[2];
        }

        xx /= edgeSubVolumesPoints.size();
        xy /= edgeSubVolumesPoints.size();
        xz /= edgeSubVolumesPoints.size();
        yy /= edgeSubVolumesPoints.size();
        yz /= edgeSubVolumesPoints.size();
        zz /= edgeSubVolumesPoints.size();

        vector weightedDir(vector::zero);

        scalar detX(yy*zz-yz*yz);
        vector axisDirX(detX,xz*yz-xy*zz,xy*yz-xz*yy);
        scalar weightX(detX*detX);
        if((weightedDir & axisDirX) < 0.0)
            weightX = -weightX;
        weightedDir += axisDirX * weightX;

        scalar detY(xx*zz-xz*xz);
        vector axisDirY(xz*yz-xy*zz,detY,xy*xz-yz*xx);
        scalar weightY(detY*detY);
        if((weightedDir & axisDirY) < 0.0)
            weightY = -weightY;
        weightedDir += axisDirY * weightY;

        scalar detZ(xx*yy-xy*xy);
        vector axisDirZ(xy*yz-xz*yy,xy*xz-yz*xx,detZ);
        scalar weightZ(detZ*detZ);
        if((weightedDir & axisDirZ) < 0.0)
            weightZ = -weightZ;
        weightedDir += axisDirZ * weightZ;

        if(mag(weightedDir) > SMALL)
        {
            normOk = true;
            normalVec = weightedDir/(mag(weightedDir)+SMALL);
        }
        if (!normOk || mag(normalVec) < 1)
            normalVec = normalVector;

        // create best fitting plane
        plane bestFitPlane(contactCenter_, normalVec);
        normalVec = bestFitPlane.normal();
        DynamicPointList commCellsPosInPlane;
        forAll (edgeSubVolumesPoints,cell)
        {
            commCellsPosInPlane.append(bestFitPlane.nearestPoint(edgeSubVolumesPoints[cell]));
		}

        vector q1(1.0, 0.0, 0.0);
        vector q2(0.0, 1.0, 0.0);
        if (abs(q1 & bestFitPlane.normal()) > abs(q2 & bestFitPlane.normal()))
            q1 = q2;

        vector u(bestFitPlane.normal() ^ q1);
        vector v(bestFitPlane.normal() ^ u);

        DynamicList<plane> clockwisePlanes;
        List<scalar> helpList(6);
        helpList[0] = 0.0;
        helpList[1] = 1.0;
        helpList[2] = 0.0;
        helpList[3] = -1.0;
        helpList[4] = 0.0;
        helpList[5] = 1.0;

        // loop over parts of plane to find and sort points
        DynamicVectorList commCellsInSections;
        for (label i = 0; i < 4; i++)
        //~ for (label i = 0; i < 8; i++)
        {
            scalar uStep(helpList[i + 1] -  helpList[i]);
            scalar vStep(helpList[i + 2] -  helpList[i + 1]);
            DynamicPointList pointsInSection;
            //~ for (scalar j = 0.0; j < 5.0; j += 1.0)
            for (scalar j = 0.0; j < 15.0; j += 1.0)
            {
                plane uPlane(contactCenter_, u*(helpList[i] + uStep*j/15.0) + v*(helpList[i + 1] + vStep*j/15.0));
                plane vPlane(contactCenter_, u*(helpList[i] + uStep*(j+1)/15.0) + v*(helpList[i + 1] + vStep*(j+1)/15.0));

                forAll (commCellsPosInPlane, celli)
                {
                    if (uPlane.sideOfPlane(commCellsPosInPlane[celli]) == 0 && vPlane.sideOfPlane(commCellsPosInPlane[celli]) == 1)
                    {
                        pointsInSection.append(commCellsPosInPlane[celli]);
                    }
                }

                if (pointsInSection.size() > SMALL)
                {
                    vector max(vector::zero);
                    scalar distance(0);
                    forAll (pointsInSection, pointI)
                    {
                         scalar magnitude = mag(pointsInSection[pointI]-contactCenter_);
                         if (magnitude > distance)
                         {
							 magnitude = distance;
							 max = pointsInSection[pointI];
						 }
					}
                    commCellsInSections.append(max);
                }
            }
        }
        // calculate contact area
        for (label i = 0; i + 1 < commCellsInSections.size(); i++)
        {
            //~ Info << commCellsInSections[i] << endl;
            vector AC(commCellsInSections[i] - contactCenter_);
            vector BC(commCellsInSections[i + 1] - contactCenter_);
            vector crossPr(AC ^ BC);
            area += mag(crossPr)/2;
        }

        if (commCellsInSections.size() > 2)
        {
            vector AC(commCellsInSections[commCellsInSections.size() - 1] - contactCenter_);
            vector BC(commCellsInSections[0] - contactCenter_);
            vector crossPr(AC ^ BC);
            area += mag(crossPr)/2;
        }
    }

    if (normalVec == vector::zero)
    {
        normalVec = normalVector;
	}

	if ((normalVector & normalVec) < 0)
	{
		normalVec = normalVec * (-1);
	}

    Tuple2<scalar,vector> returnValue(area,normalVec);

    return returnValue;
}
//---------------------------------------------------------------------------//
std::vector<sub_contact> virtual_mesh::findsub_contacts(sub_volume& sV)
{
    if (sV.hasChildSubVolumes())
    {
        std::vector<sub_contact> childsSubContacts;
        for (sub_volume& child : sV.childSubVolumes())
        {
            std::vector<sub_contact> childSubContacts(findsub_contacts(child));
            std::move(childSubContacts.begin(), childSubContacts.end(), std::back_inserter(childsSubContacts));
        }

        if (childsSubContacts.size() == 0)
        {
            return {};
        }

        std::vector<sub_contact> returnValue;
        returnValue.push_back(childsSubContacts.back());
        childsSubContacts.pop_back();

        while (childsSubContacts.size() > 0)
        {
            bool found = false;
            auto childSCIter = childsSubContacts.begin();
            while (childSCIter != childsSubContacts.end())
            {
                if (canCombineSubContacts(
                    returnValue.back(),
                    *childSCIter
                ))
                {
                    sub_contact& childSC = *childSCIter;
                    for (auto& sV : childSC.getSubVolumes())
                    {
                        returnValue.back().addSubVolume(sV);
                    }

                    childsSubContacts.erase(childSCIter);
                    found = true;
                    break;
                }

                ++childSCIter;
            }

            if (!found)
            {
                returnValue.push_back(childsSubContacts.back());
                childsSubContacts.pop_back();
            }
        }

        return returnValue;
    }

    std::vector<sub_contact> returnValue;
    if (sV.cVolumeInfo().volumeType_ != volumeType::outside
        && sV.tVolumeInfo().volumeType_ != volumeType::outside)
    {
        returnValue.push_back(sub_contact());
        returnValue.back().addSubVolume(std::make_shared<sub_volume>(sV));
    }

    return returnValue;
}
//---------------------------------------------------------------------------//
bool virtual_mesh::canCombineSubContacts(
    sub_contact& main,
    sub_contact& comp
) const
{
    if (!main.getBounds().overlaps(comp.getBounds()))
    {
        return false;
    }

    for (const auto& sV : comp.getSubVolumes())
    {
        if (main.canCombine(*sV))
        {
            return true;
        }
    }

    return false;
}
//---------------------------------------------------------------------------//
// ************************************************************************* //
