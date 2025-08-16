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
#include "line_interpolInfo.H"

#define ORDER 2

using namespace Foam;

//---------------------------------------------------------------------------//
line_interpolInfo::line_interpolInfo
(
    const  fvMesh&   mesh,
    std::shared_ptr<geom_model>& gModel
)
:
interpol_info(mesh, gModel)
{}
line_interpolInfo::~line_interpolInfo()
{}
//---------------------------------------------------------------------------//
void line_interpolInfo::setIntpInfo()
{
    const DynamicLabelList& cSurfCells = getSurfCells();

    resetIntpInfo(cSurfCells.size());
    List<point>& ibPoints = getIbPoints();
    List<vector>& ibNormals = getIbNormals();
    List<List<interpol_point>>& interpol_points = getIntPoints();


    forAll (cSurfCells, cellI)
    {
        label scell = cSurfCells[cellI];
        scalar intDist = Foam::pow(mesh_.V()[scell],0.333);

        geom_model_->getClosestPointAndNormal(
            mesh_.C()[scell],
            intDist*2*vector::one,
            ibPoints[cellI],
            ibNormals[cellI]
        );

        interpol_points[cellI].setSize(ORDER);
        interpol_point cIntPoint
        (
            ibPoints[cellI],
            scell,
            Pstream::myProcNo()
        );
        point cPoint;

        for(label i = 0; i < ORDER; ++i)
        {
            cPoint = cIntPoint.iPoint_;
            do {
                cPoint += ibNormals[cellI]*intDist;
            } while(pointInCell(cPoint, cIntPoint.iCell_));

            interpol_points[cellI][i] = findIntPoint(cIntPoint, cPoint);
            correctIntPoint(ibPoints[cellI], interpol_points[cellI][i]);
            cIntPoint = interpol_points[cellI][i];

            if(cIntPoint.iProc_ != Pstream::myProcNo())
            {
                break;
            }
        }
    }
    syncIntPoints();
}
//---------------------------------------------------------------------------//
void line_interpolInfo::correctIntPoint
(
    point ibPoint,
    interpol_point& cPoint
)
{
    if(cPoint.iProc_ != Pstream::myProcNo())
    {
        return;
    }

    vector closestPoint = getClosestPoint(ibPoint, cPoint);

    if(pointInCell(closestPoint, cPoint.iCell_))
    {
        cPoint.iPoint_ = closestPoint;
    }
    else
    {
        const labelList& cellFaces(mesh_.cells()[cPoint.iCell_]);

        forAll (cellFaces, fi)
        {
            const face faceI = mesh_.faces()[cellFaces[fi]];
            vector dir = closestPoint - cPoint.iPoint_;

            pointHit pHit = faceI.ray(
                cPoint.iPoint_,
                dir,
                mesh_.points()
            );

            if(pHit.hit())
            {
                vector newP = 0.95*(pHit.hitPoint() - cPoint.iPoint_);
                newP += cPoint.iPoint_;

                if(pointInCell(newP, cPoint.iCell_))
                {
                    cPoint.iPoint_ = newP;
                    break;
                }
            }
        }
    }
}
//---------------------------------------------------------------------------//
vector line_interpolInfo::getClosestPoint
(
    vector ibPoint,
    interpol_point& cPoint
)
{
    vector dir = cPoint.iPoint_ - ibPoint;
    dir /= mag(dir);

    vector dirToC = mesh_.C()[cPoint.iCell_] - ibPoint;

    return ibPoint + dir*(dirToC&dir);
}
//---------------------------------------------------------------------------//
interpol_point line_interpolInfo::findIntPoint
(
    interpol_point& fromP,
    point& endP
)
{
    if(fromP.iPoint_ == endP)
    {
        return interpol_point();
    }

    interpol_point retP
    (
        endP,
        fromP.iCell_,
        fromP.iProc_
    );

    if(fromP.iProc_ == Pstream::myProcNo())
    {
        label faceInDir = -1;
        while(!pointInCell(retP.iPoint_, retP.iCell_))
        {
            faceInDir = getFaceInDir(retP, faceInDir);
            if (!mesh_.isInternalFace(faceInDir))
            {
                label facePatchId(mesh_.boundaryMesh().whichPatch(faceInDir));
                const polyPatch& cPatch = mesh_.boundaryMesh()[facePatchId];

                if (cPatch.type() == "processor")
                {
                    const processorPolyPatch& procPatch
                        = refCast<const processorPolyPatch>(cPatch);

                    label sProc = (Pstream::myProcNo() == procPatch.myProcNo())
                        ? procPatch.neighbProcNo() : procPatch.myProcNo();

                    retP.iCell_ = cPatch.whichFace(faceInDir);
                    retP.iProc_ = sProc;

                    return retP;
                }
                else
                {
                    retP.iProc_ = -1;
                    return retP;
                }
            }

            label owner(mesh_.owner()[faceInDir]);
            label neighbour(mesh_.neighbour()[faceInDir]);
            retP.iCell_ = (retP.iCell_ == neighbour) ? owner : neighbour;
        }

        return retP;
    }

    return retP;
}
//---------------------------------------------------------------------------//
label line_interpolInfo::getFaceInDir
(
    const interpol_point& retPoint,
    const label prevFace
)
{
    label faceToReturn = -1;
    vector dir = retPoint.iPoint_ - mesh_.C()[retPoint.iCell_];

    const labelList& cellFaces(mesh_.cells()[retPoint.iCell_]);
    scalar dotProd(-GREAT);

    forAll (cellFaces, faceI)
    {
        label fI = cellFaces[faceI];
        if(fI != prevFace)
        {
            vector outNorm = (mesh_.faceOwner()[fI] == retPoint.iCell_)
                ? mesh_.Sf()[fI] : (-1*mesh_.Sf()[fI]);

            scalar auxDotProd(outNorm & dir);
            if (auxDotProd > dotProd)
            {
                dotProd = auxDotProd;
                faceToReturn = fI;
            }
        }
    }

    return faceToReturn;
}
//---------------------------------------------------------------------------//
bool line_interpolInfo::pointInCell
(
    point pToCheck,
    label cToCheck
)
{
    const labelList& cellFaces(mesh_.cells()[cToCheck]);
    forAll(cellFaces, faceI)
    {
        label fI = cellFaces[faceI];
        vector outNorm = mesh_.Sf()[fI];
        outNorm = (mesh_.faceOwner()[fI] == cToCheck) ? outNorm : (-1*outNorm);

        if (((pToCheck - mesh_.Cf()[fI]) & outNorm) > 0)
        {
            return false;
        }
    }
    return true;
}
//---------------------------------------------------------------------------//
void line_interpolInfo::syncIntPoints()
{
    List<point>& ibPoints = getIbPoints();
    List<List<interpol_point>>& interpol_points = getIntPoints();

    List<DynamicPointList> ibPointsToSync(Pstream::nProcs());
    List<DynamicPointList> interpol_pointToSync(Pstream::nProcs());
    List<DynamicLabelList> faceLabelToSync(Pstream::nProcs());
    List<DynamicLabelList> orderToSync(Pstream::nProcs());
    List<DynamicLabelList> labelToSync(Pstream::nProcs());

    forAll(ibPoints, pI)
    {
        forAll(interpol_points[pI], ipI)
        {
            if(interpol_points[pI][ipI].iProc_ != Pstream::myProcNo()
                &&
                interpol_points[pI][ipI].iProc_ != -1)
            {
                interpol_point& cIntPoint = interpol_points[pI][ipI];
                ibPointsToSync[cIntPoint.iProc_].append(ibPoints[pI]);
                interpol_pointToSync[cIntPoint.iProc_].append(cIntPoint.iPoint_);
                faceLabelToSync[cIntPoint.iProc_].append(cIntPoint.iCell_);
                orderToSync[cIntPoint.iProc_].append(ipI);
                labelToSync[cIntPoint.iProc_].append(pI);
            }
        }
    }

    PstreamBuffers pBufsIbP(Pstream::commsTypes::nonBlocking);
    PstreamBuffers pBufsIntP(Pstream::commsTypes::nonBlocking);
    PstreamBuffers pBufsFaceL(Pstream::commsTypes::nonBlocking);
    PstreamBuffers pBufsOrder(Pstream::commsTypes::nonBlocking);
    PstreamBuffers pBufsLabel(Pstream::commsTypes::nonBlocking);

    for (label proci = 0; proci < Pstream::nProcs(); proci++)
    {
        if(proci != Pstream::myProcNo())
        {
            UOPstream sendIbP(proci, pBufsIbP);
            UOPstream sendIntP(proci, pBufsIntP);
            UOPstream sendFaceL(proci, pBufsFaceL);
            UOPstream sendOrder(proci, pBufsOrder);
            UOPstream sendLabel(proci, pBufsLabel);

            sendIbP << ibPointsToSync[proci];
            sendIntP << interpol_pointToSync[proci];
            sendFaceL << faceLabelToSync[proci];
            sendOrder << orderToSync[proci];
            sendLabel << labelToSync[proci];
        }
    }

    pBufsIbP.finishedSends();
    pBufsIntP.finishedSends();
    pBufsFaceL.finishedSends();
    pBufsOrder.finishedSends();
    pBufsLabel.finishedSends();

    List<DynamicPointList> ibPointsRecv(Pstream::nProcs());
    List<DynamicPointList> interpol_pointRecv(Pstream::nProcs());
    List<DynamicLabelList> faceLabelRecv(Pstream::nProcs());
    List<DynamicLabelList> orderRecv(Pstream::nProcs());
    List<DynamicLabelList> labelRecv(Pstream::nProcs());

    for (label proci = 0; proci < Pstream::nProcs(); proci++)
    {
        if (proci != Pstream::myProcNo())
        {
            UIPstream recvIbP(proci, pBufsIbP);
            UIPstream recvIntP(proci, pBufsIntP);
            UIPstream recvFaceL(proci, pBufsFaceL);
            UIPstream recvOrder(proci, pBufsOrder);
            UIPstream recvLabel(proci, pBufsLabel);

            DynamicPointList recIbP (recvIbP);
            DynamicPointList recIntP (recvIntP);
            DynamicLabelList recFaceL (recvFaceL);
            DynamicLabelList recOrder (recvOrder);
            DynamicLabelList recLabel (recvLabel);

            ibPointsRecv[proci] = recIbP;
            interpol_pointRecv[proci] = recIntP;
            faceLabelRecv[proci] = recFaceL;
            orderRecv[proci] = recOrder;
            labelRecv[proci] = recLabel;
        }
    }

    pBufsIbP.clear();
    pBufsIntP.clear();
    pBufsFaceL.clear();
    pBufsOrder.clear();
    pBufsLabel.clear();

    List<DynamicLabelList> cellLabelRecv(Pstream::nProcs());

    forAll (mesh_.boundaryMesh(), patchi)
    {
        const polyPatch& cPatch = mesh_.boundaryMesh()[patchi];
        if (cPatch.type() == "processor")
        {
            const processorPolyPatch& procPatch
                = refCast<const processorPolyPatch>(cPatch);

            label sProc = (Pstream::myProcNo() == procPatch.myProcNo())
                ? procPatch.neighbProcNo() : procPatch.myProcNo();

            cellLabelRecv[sProc].setSize(faceLabelRecv[sProc].size());
            forAll(faceLabelRecv[sProc], faceI)
            {
                cellLabelRecv[sProc][faceI]
                    = mesh_.faceOwner()[cPatch.start()
                    + faceLabelRecv[sProc][faceI]];
            }
        }
    }

    List<DynamicPointList> interpol_pointToRetr(Pstream::nProcs());
    List<DynamicLabelList> intCellToRetr(Pstream::nProcs());
    List<DynamicLabelList> intProcToRetr(Pstream::nProcs());
    List<DynamicLabelList> orderToRetr(Pstream::nProcs());
    List<DynamicLabelList> labelToRetr(Pstream::nProcs());

    forAll(ibPointsRecv, proci)
    {
        forAll(ibPointsRecv[proci], ibpI)
        {
            point cPoint = mesh_.C()[cellLabelRecv[proci][ibpI]];
            interpol_point cIntPoint
            (
                cPoint,
                cellLabelRecv[proci][ibpI],
                proci
            );

            interpol_point foundP =
                findIntPoint(cIntPoint, interpol_pointRecv[proci][ibpI]);

            scalar intDist = Foam::pow(mesh_.V()[foundP.iCell_],0.333);
            vector dir = foundP.iPoint_ - ibPointsRecv[proci][ibpI];
            dir /= mag(dir);

            correctIntPoint(ibPointsRecv[proci][ibpI], foundP);

            interpol_pointToRetr[proci].append(foundP.iPoint_);
            intCellToRetr[proci].append(foundP.iCell_);
            intProcToRetr[proci].append(foundP.iProc_);
            orderToRetr[proci].append(orderRecv[proci][ibpI]);
            labelToRetr[proci].append(labelRecv[proci][ibpI]);

            cIntPoint = foundP;

            for(label i = orderRecv[proci][ibpI] + 1; i < ORDER; ++i)
            {
                cPoint = cIntPoint.iPoint_;
                do {
                    cPoint += dir*intDist;
                } while(pointInCell(cPoint, cIntPoint.iCell_));

                foundP = findIntPoint(cIntPoint, cPoint);
                correctIntPoint(ibPointsRecv[proci][ibpI], foundP);
                cIntPoint = foundP;

                if(cIntPoint.iProc_ != proci)
                {
                    break;
                }

                interpol_pointToRetr[proci].append(foundP.iPoint_);
                intCellToRetr[proci].append(foundP.iCell_);
                intProcToRetr[proci].append(foundP.iProc_);
                orderToRetr[proci].append(i);
                labelToRetr[proci].append(labelRecv[proci][ibpI]);
            }
        }
    }

    for (label proci = 0; proci < Pstream::nProcs(); proci++)
    {
        if(proci != Pstream::myProcNo())
        {
            UOPstream sendIbP(proci, pBufsIbP);
            UOPstream sendIntP(proci, pBufsIntP);
            UOPstream sendFaceL(proci, pBufsFaceL);
            UOPstream sendOrder(proci, pBufsOrder);
            UOPstream sendLabel(proci, pBufsLabel);

            sendIbP << intProcToRetr[proci];
            sendIntP << interpol_pointToRetr[proci];
            sendFaceL << intCellToRetr[proci];
            sendOrder << orderToRetr[proci];
            sendLabel << labelToRetr[proci];
        }
    }

    pBufsIbP.finishedSends();
    pBufsIntP.finishedSends();
    pBufsFaceL.finishedSends();
    pBufsOrder.finishedSends();
    pBufsLabel.finishedSends();

    List<DynamicPointList> interpol_pointCmpl(Pstream::nProcs());
    List<DynamicLabelList> intCellCmpl(Pstream::nProcs());
    List<DynamicLabelList> intProcCmpl(Pstream::nProcs());
    List<DynamicLabelList> orderCmpl(Pstream::nProcs());
    List<DynamicLabelList> labelCmpl(Pstream::nProcs());

    for (label proci = 0; proci < Pstream::nProcs(); proci++)
    {
        if (proci != Pstream::myProcNo())
        {
            UIPstream recvIntP(proci, pBufsIntP);
            UIPstream recvCell(proci, pBufsFaceL);
            UIPstream recvProc(proci, pBufsIbP);
            UIPstream recvOrder(proci, pBufsOrder);
            UIPstream recvLabel(proci, pBufsLabel);

            DynamicPointList recIntP (recvIntP);
            DynamicLabelList recCell (recvCell);
            DynamicLabelList recProc (recvProc);
            DynamicLabelList recOrder (recvOrder);
            DynamicLabelList recLabel (recvLabel);

            interpol_pointCmpl[proci] = recIntP;
            intCellCmpl[proci] = recCell;
            intProcCmpl[proci] = recProc;
            orderCmpl[proci] = recOrder;
            labelCmpl[proci] = recLabel;
        }
    }

    pBufsIbP.clear();
    pBufsIntP.clear();
    pBufsFaceL.clear();
    pBufsOrder.clear();
    pBufsLabel.clear();

    forAll(interpol_pointCmpl, proci)
    {
        forAll(interpol_pointCmpl[proci], iPointI)
        {
            interpol_point cIntPoint
            (
                interpol_pointCmpl[proci][iPointI],
                intCellCmpl[proci][iPointI],
                intProcCmpl[proci][iPointI]
            );

            interpol_points[labelCmpl[proci][iPointI]][orderCmpl[proci][iPointI]]
                = cIntPoint;
        }
    }
}
//---------------------------------------------------------------------------//
