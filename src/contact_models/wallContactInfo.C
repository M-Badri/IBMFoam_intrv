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
#include "wallContactInfo.H"

#include "inter_adhesion.H"
#include "wall_mat_info.H"

#include "virtual_meshLevel.H"
#include "contact_model_info.H"

using namespace Foam;
//---------------------------------------------------------------------------//
wallContactInfo::wallContactInfo
(
    ibContactClass& cClass,
    ibContactVars& cVars
)
:
ibContactClass_(cClass),
ibContactVars_(cVars)
{
    body_id_ = cVars.body_id_;
    const material_info& cMatInfo(cClass.getMatInfo());

    const List<string> cntPatches = wall_mat_info::getWallPatches();
    forAll(cntPatches, patchI)
    {
        const material_info& wInfo = wall_mat_info::getWallMatInfo()[cntPatches[patchI]];

        scalar adhPot = 0;
        string adhPotKey;

        if(cMatInfo.getMaterial() < wInfo.getMaterial())
        {
            adhPotKey = cMatInfo.getMaterial();
            adhPotKey += "-";
            adhPotKey += wInfo.getMaterial();
        }
        else
        {
            adhPotKey = wInfo.getMaterial();
            adhPotKey += "-";
            adhPotKey += cMatInfo.getMaterial();
        }
        if(inter_adhesion::getInterAdhesion().found(adhPotKey))
        {
            adhPot = inter_adhesion::getInterAdhesion()[adhPotKey];
        }

        scalar aY = 1/((1 - sqr(cMatInfo.getNu()))/cMatInfo.getY()
            + (1 - sqr(wInfo.getNu()))/wInfo.getY());
        scalar aG = 1/(2*(2 - cMatInfo.getNu())*(1 + cMatInfo.getNu())
            /cMatInfo.getY() + 2*(2 - wInfo.getNu())
            *(1 + wInfo.getNu())/wInfo.getY());
        scalar aMu = (cMatInfo.getMu()+wInfo.getMu())/2;
        scalar maxAdhN = cMatInfo.getAdhN() + wInfo.getAdhN() - 2*adhPot;
        scalar reduceBeta =
        (
           (-1)*sqrt(5.0)*log((0.5*(cMatInfo.getEps()+wInfo.getEps())))/
            (sqrt(sqr(log((0.5*(cMatInfo.getEps()+wInfo.getEps()))))+
            sqr(Foam::constant::mathematical::pi)))
        );
        wallMeanPars_.insert(
            cntPatches[patchI],
            physicalProperties(aY, aG, aMu, maxAdhN, 0, 0, reduceBeta)
        );

    }

    reduceM_ = 0;
}

wallContactInfo::~wallContactInfo()
{
}
//---------------------------------------------------------------------------//
boundBox wallContactInfo::constructBoundBox
(
    boundBox& bodyBB
)
{
    vector minPoint = bodyBB.min();
    vector maxPoint = bodyBB.max();

    for(int i=0;i<3;i++)
    {
        minPoint[i] = floor(minPoint[i]/virtual_meshLevel::getCharCellSize())*virtual_meshLevel::getCharCellSize();
        maxPoint[i] = ceil(maxPoint[i]/virtual_meshLevel::getCharCellSize())*virtual_meshLevel::getCharCellSize();
    }
    boundBox BB = boundBox(minPoint,maxPoint);
    return BB;
}
//---------------------------------------------------------------------------//
void wallContactInfo::constructSM()
{
    boundBox bodyBB = ibContactClass_.get_geom_model().getBounds();

    boundBox BB = constructBoundBox(bodyBB);

    vector cellNVector = vector(
        ceil((BB.span()[0]/virtual_meshLevel::getCharCellSize())),
        ceil((BB.span()[1]/virtual_meshLevel::getCharCellSize())),
        ceil((BB.span()[2]/virtual_meshLevel::getCharCellSize()))
    );

    const scalar charCellSize(virtual_meshLevel::getCharCellSize());

    SM_.reset(new spectator_mesh(cellNVector,BB,charCellSize));
}
//---------------------------------------------------------------------------//
bool wallContactInfo::isInsidePlane(
    vector checkedPoint,
    const string& wall,
    const HashTable<List<vector>,string,Hash<string>>& wallPatches
)
{
    const vector& normalVector = wallPatches[wall][0];
    const vector& centerPoint = wallPatches[wall][1];
    vector testVector = checkedPoint - centerPoint;
    if((testVector & normalVector) < 0)
    {
        return true;
    }
    return false;
}
//---------------------------------------------------------------------------//
bool wallContactInfo::detectWallContact(
    const HashTable<List<vector>,string,Hash<string>>& wallPatches
)
{
    contactPatches_.clear();
    clearOldContact();
    boundBox bodyBB = ibContactClass_.get_geom_model().getBounds();
    pointField bBpoints = bodyBB.points();
    List<bool> isPatchInContact;
    const stringList wallPatchNames = wallPatches.toc();
    isPatchInContact.setSize(wallPatchNames.size(),false);
    bool bBContact(false);
    forAll (bBpoints,bBP)
    {
        forAll(wallPatchNames,wP)
        {
            if(!isInsidePlane(bBpoints[bBP], wallPatchNames[wP], wallPatches))
            {
                isPatchInContact[wP] = true;
                bBContact = true;
            }
        }
    }

    forAll(isPatchInContact,wP)
    {
        if(isPatchInContact[wP])
        {
            contactPatches_.append(wallPatchNames[wP]);
        }
    }

    return bBContact;

}
//---------------------------------------------------------------------------//
void wallContactInfo::findContactAreas()
{
    autoPtr<DynamicVectorList> contactSTLPoints(
        new DynamicVectorList);

    pointField bodyPoints = ibContactClass_.get_geom_model().getBodyPoints();
    forAll(bodyPoints,bP)
    {
        forAll(contactPatches_,wP)
        {
            if(!isInsidePlane(bodyPoints[bP],contactPatches_[wP]))
            {
                contactSTLPoints().append(bodyPoints[bP]);
            }

        }
    }


    if(contactSTLPoints().size() > SMALL)
    {
        constructSM();
        List<DynamicList<vector>> possibleSMContact = detectPossibleSMContact(contactSTLPoints(),contactPatches_);
        forAll(possibleSMContact,SC)
        {
            label emptyCells(0);
            pointField overallContactPoints;
            List<string> contactPatches;

            List<Tuple2<point,boundBox>> sMExportList;
            List<Tuple2<point,boundBox>> sMPlaneList;
            List<Tuple2<point,boundBox>> sMInternal;

            boundBox cBbox;
            forAll(possibleSMContact[SC],item)
            {
                Tuple2<point,boundBox> sMExport;
                forAll(contactPatches_,wP)
                {
                    boundBox sMEBB = correctSMBBforWall(SM_().getElementBB(possibleSMContact[SC][item]),contactPatches_[wP]);
                    if(sMEBB.minDim() != 0)
                    {
                        sMExport.first() = SM_()[possibleSMContact[SC][item]].initPoint;
                        sMExport.second() = sMEBB;

                        overallContactPoints.append(sMEBB.points());

                        if(!SM_()[possibleSMContact[SC][item]].isInternal)
                        {
                            sMExportList.append(sMExport);
                        }
                        else
                        {
                            sMInternal.append(sMExport);
                        }
                        break;
                    }
                    else
                    {
                        emptyCells++;
                    }
                }
            }
            forAll(contactPatches_,wP)
            {

                boundBox sCBBox = getSCBBox(possibleSMContact[SC]);
                List<bool> isPathechInSC;
                isPathechInSC.setSize(contactPatches_.size(),false);
                pointField sCBBPoints = sCBBox.points();
                forAll(sCBBPoints,bBP)
                {
                    if(!isInsidePlane(sCBBPoints[bBP],contactPatches_[wP]))
                    {
                        isPathechInSC[wP]+= true;
                    }
                }
                // InfoH << DEM_Info << "-- SM sCBBox "<< sCBBox << endl;
                if(isPathechInSC[wP])
                {

                    cBbox = boundBox(overallContactPoints,false);
                    boundBox planeBox = contactPlaneBBox(constructVMBox(cBbox,contactPatches_[wP]),contactPatches_[wP]);

                    contactPatches.append(contactPatches_[wP]);

                    forAll(contactPatches_,wP2)
                    {
                        if(wP != wP2)
                        {
                            forAll(sCBBPoints,bBP)
                            {
                                if(!isInsidePlane(sCBBPoints[bBP],contactPatches_[wP2]))
                                {
                                    isPathechInSC[wP2]+= true;
                                }
                            }

                            if(isPathechInSC[wP2])
                            {
                                pointField bbPoints = planeBox.points();
 
                                pointField newBBPoints;
                                forAll(bbPoints,bBP)
                                {
                                    if(!isInsidePlane(bbPoints[bBP],contactPatches_[wP2]))
                                    {
                                        newBBPoints.append(getPlanePoint(bbPoints[bBP],contactPatches_[wP2]));
                                    }
                                    else
                                    {
                                        newBBPoints.append(bbPoints[bBP]);
                                    }
                                }
                                planeBox = boundBox(newBBPoints,false);
                            }

                        }
                    }
 
                    Tuple2<point,boundBox> sMPlane(planeBox.midpoint(),planeBox);
                    sMPlaneList.append(sMPlane);
                }

            }

            setNewSubContact(sMExportList,sMPlaneList,contactPatches,sMInternal,cBbox);
        }
    }
}
//---------------------------------------------------------------------------//
List<DynamicList<vector>> wallContactInfo::detectPossibleSMContact
(
    DynamicVectorList& contactPoints,
    List<string>& contactPatches
)
{
    vectorHashSet possibleContactElements;
    vectorHashSet checkedOctreeFaces;
    List<DynamicVectorList> baseSubContactList;
    label iterMax(SM_().matrixSize_[0]*SM_().matrixSize_[1]*SM_().matrixSize_[2]);

    autoPtr<DynamicVectorList> contactElements(
        new DynamicVectorList);



    forAll(contactPoints,cP)
    {
        vector elementIndex(SM_().getSMCentroidIndex(contactPoints[cP]));
        if(!possibleContactElements.found(elementIndex))
        {
            possibleContactElements.insert(elementIndex);
            SM_()[elementIndex].stlFound = true;
            SM_()[elementIndex].initPointSet = true;
            SM_()[elementIndex].initPoint = contactPoints[cP];
            contactElements().append(elementIndex);
        }
    }

    autoPtr<DynamicVectorList> nextToCheck(
        new DynamicVectorList);

    autoPtr<DynamicVectorList> auxToCheck(
            new DynamicVectorList);

    autoPtr<DynamicVectorList> sub_contactElements(
            new DynamicVectorList);

    label iterCount2(0);
   
    while((contactElements().size() > SMALL) and iterCount2++ < iterMax)
    {
        label iterCount(0);
        nextToCheck().clear();
        sub_contactElements().clear();

        nextToCheck().append(contactElements()[0]);

        while ((nextToCheck().size() > 0) and iterCount++ < iterMax)
        {
            auxToCheck().clear();
            forAll(nextToCheck(),cE)
            {
                if (SM_()[nextToCheck()[cE]].toCheck)
                {
                    checkSMElement(nextToCheck()[cE],contactPatches);
                    bool isInMesh(false);
                    bool isInBody(false);
                    bool allInMesh(false);
                    bool allInBody(false);
                    checkElement(nextToCheck()[cE],isInMesh,isInBody,allInMesh,allInBody);
                    if(isInBody && isInMesh)
                    {
                        if(!allInMesh)
                        {
                            sub_contactElements().append(nextToCheck()[cE]);
                            auxToCheck().append(SM_().faceNeighbourElements(nextToCheck()[cE]));
                            checkedOctreeFaces.insert(nextToCheck()[cE]);
                            if(allInBody)
                            {
                                SM_()[nextToCheck()[cE]].isInternal = allInBody;
                            }
                        }
                    }
                    if((isInBody && !isInMesh)|| SM_()[nextToCheck()[cE]].stlFound)
                    {
                        sub_contactElements().append(nextToCheck()[cE]);
                        checkedOctreeFaces.insert(nextToCheck()[cE]);
                        auxToCheck().append(SM_().faceNeighbourElements(nextToCheck()[cE]));
                        if(allInBody)
                        {
                            SM_()[nextToCheck()[cE]].isInternal = allInBody;
                        }
                    }
                }
            }
            const autoPtr<DynamicVectorList> helpPtr(nextToCheck.ptr());
            nextToCheck.set(auxToCheck.ptr());
            auxToCheck = helpPtr;
        }
        if(sub_contactElements().size() > 0)
        {
            baseSubContactList.append(sub_contactElements());
        }
        for (auto it = contactElements().begin(); it != contactElements().end();)
        {
            if (checkedOctreeFaces.found(*it))
            {
                it = contactElements().erase(it);
            }
            else
            {
                ++it;
            }
        }
    }
    return baseSubContactList;
}
//---------------------------------------------------------------------------//
boundBox wallContactInfo::getSCBBox
(
    DynamicVectorList& sub_contactAreas
)
{
    pointField contactArea;
    forAll(sub_contactAreas,sCE)
    {
        contactArea.append(SM_()[sub_contactAreas[sCE]].center);
        List<vector> vertexList = SM_().elementVertexIndexies(sub_contactAreas[sCE]);
        forAll(vertexList,vL)
        {
            contactArea.append(SM_()(vertexList[vL]).center);
        }
    }
    return boundBox(contactArea,false);
}
//---------------------------------------------------------------------------//
boundBox wallContactInfo::constructVMBox
(
    boundBox& baseContactAreaBB,
    string& wallName
)
{

    pointField bbPoints = baseContactAreaBB.points();
    pointField newBBPoints;
    forAll(bbPoints,bP)
    {
        if(isInsidePlane(bbPoints[bP],wallName))
        {
            newBBPoints.append(getPlanePoint(bbPoints[bP],wallName));
        }
        else
        {
           newBBPoints.append(bbPoints[bP]);
        }
    }
    return boundBox(newBBPoints,false);
}
//---------------------------------------------------------------------------//
vector wallContactInfo::getPlanePoint
(
    vector pointInDomain,
    string wallName
)
{

    const vector& normalVector = wall_plane_info::getWallPlaneInfo()[wallName][0];
    const vector& centerPoint = wall_plane_info::getWallPlaneInfo()[wallName][1];


    scalar fracI = normalVector & (centerPoint - pointInDomain);
    scalar fracII = magSqr(normalVector);


    return pointInDomain + (fracI/fracII)*normalVector;
}
//---------------------------------------------------------------------------//
void wallContactInfo::checkSMElement
(
    vector& index,
    List<string> contactPatches
)
{
    label iC(0);
    if(SM_()[index].toCheck)
    {
        point centroidPoint = SM_()[index].center;
        SM_()[index].toCheck = false;
        SM_()[index].isCBody = ibContactClass_.get_geom_model().pointInside(centroidPoint);
        bool isMeshLocC(true);
        while(iC < contactPatches.size() and !SM_()[index].isMesh)
        {

            isMeshLocC *= isInsidePlane(centroidPoint,contactPatches[iC]);
            iC++;
        }

        SM_()[index].isMesh = isMeshLocC;

        const List<vector> verticesLabels = SM_().elementVertexIndexies(index);
        forAll(verticesLabels, vL)
        {
            iC = 0;
            if(SM_()(verticesLabels[vL]).toCheck)
            {
                point vertexPoint = SM_()(verticesLabels[vL]).center;
                SM_()(verticesLabels[vL]).toCheck = false;
                SM_()(verticesLabels[vL]).isCBody = ibContactClass_.get_geom_model().pointInside(vertexPoint);
                bool isMeshLocV(true);
                while(iC < contactPatches.size() and !SM_()(verticesLabels[vL]).isMesh)
                {
                    isMeshLocV *= isInsidePlane(vertexPoint,contactPatches[iC]);
                    iC ++;
                }
                SM_()(verticesLabels[vL]).isMesh = isMeshLocV;
            }
        }
    }
}
//---------------------------------------------------------------------------//
void wallContactInfo::checkElement
(
    vector& index,
    bool& inMesh,
    bool& inBody,
    bool& allInMesh,
    bool& allInBody
)
{
    bool inMeshLocal(false);
    bool inBodyLocal(false);
    bool allInMeshLocal(true);
    bool allInBodyLocal(true);

    inMeshLocal =  inMeshLocal || SM_()[index].isMesh;
    allInMeshLocal *= SM_()[index].isMesh;
    inBodyLocal =  inBodyLocal || SM_()[index].isCBody;
    allInBodyLocal *= SM_()[index].isCBody;
    if(SM_()[index].isCBody && !SM_()[index].initPointSet)
    {
        SM_()[index].initPoint = SM_()[index].center;
        SM_()[index].initPointSet = true;
    }
    const List<vector> verticesLabels = SM_().elementVertexIndexies(index);

    forAll(verticesLabels, vL)
    {
        inMeshLocal =  inMeshLocal || SM_()(verticesLabels[vL]).isMesh;
        allInMeshLocal *= SM_()(verticesLabels[vL]).isMesh;
        inBodyLocal =  inBodyLocal || SM_()(verticesLabels[vL]).isCBody;
        allInBodyLocal *= SM_()(verticesLabels[vL]).isCBody;
        if(SM_()(verticesLabels[vL]).isCBody && !SM_()[index].isCBody && !SM_()[index].initPointSet)
        {
            SM_()[index].initPoint = SM_()(verticesLabels[vL]).center;
            SM_()[index].initPointSet = true;
        }
    }

    inMesh = inMeshLocal;
    inBody = inBodyLocal;
    allInMesh = allInMeshLocal;
    allInBody = allInBodyLocal;
}
//---------------------------------------------------------------------------//
boundBox wallContactInfo::contactPlaneBBox
(
    boundBox contactBoundBox,
    string wallName
)
{
    pointField planeBB;
    const pointField contactBB = contactBoundBox.points();
    forAll(contactBB,bBP)
    {
        planeBB.append(getPlanePoint(contactBB[bBP],wallName));
    }
    return boundBox(planeBB,false);
}
//---------------------------------------------------------------------------//
boundBox wallContactInfo::correctSMBBforWall
(
    boundBox bB,
    string wallName
)
{
    pointField newBB;
    const pointField elementBBPoints = bB.points();
    forAll(elementBBPoints,bBP)
    {
        if(isInsidePlane(elementBBPoints[bBP],wallName))
        {
            newBB.append(getPlanePoint(elementBBPoints[bBP],wallName));
        }
        else
        {
           newBB.append(elementBBPoints[bBP]);
        }
    }
    return boundBox(newBB,false);
}
//---------------------------------------------------------------------------//
void wallContactInfo::clearOldContact()
{
    subCList_.clear();
}
//---------------------------------------------------------------------------//
void wallContactInfo::setNewSubContact(
    const List<Tuple2<point,boundBox>>& contactBBData,
    const List<Tuple2<point,boundBox>>& planeBBData,
    const List<string>& contactPatches,
    const List<Tuple2<point,boundBox>>& isInternal,
    boundBox BB
)
{
    subCList_.emplace_back(std::make_shared<wallSubContactInfo>
    (contactBBData, planeBBData,contactPatches,isInternal,wallMeanPars_,BB,body_id_));

}
//---------------------------------------------------------------------------//
void wallContactInfo::registerSubContactList(DynamicList<wallSubContactInfo*>& wall_contact_list)
{
    for(auto sC : subCList_)
    {
        wall_contact_list.append(sC.get());
    }
}
//---------------------------------------------------------------------------//
