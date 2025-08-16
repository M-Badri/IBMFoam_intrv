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
#include "line_interpol.H"

#define ORDER 2

using namespace Foam;

//---------------------------------------------------------------------------//
line_interpol::line_interpol(dictionary& interpDict):
interpDict_(interpDict)
{}
line_interpol::~line_interpol()
{}
//---------------------------------------------------------------------------//
void line_interpol::ib_interpolate
(
    interpol_info& intpInfo,
    volVectorField& Ui,
    vectorField ibPointsVal,
    const Foam::fvMesh& mesh
)
{
    line_interpolInfo& lsInfo
        = dynamic_cast<line_interpolInfo&>(intpInfo);

    correctVelocity
    (
        lsInfo,
        Ui,
        ibPointsVal,
        mesh
    );
}
//---------------------------------------------------------------------------//
void line_interpol::correctVelocity
(
    line_interpolInfo& intpInfo,
    volVectorField& Ui,
    vectorField& ibPointsVal,
    const Foam::fvMesh& mesh
)
{
    const DynamicLabelList& cSurfCells = intpInfo.getSurfCells();

    List<point>& ibPoints = intpInfo.getIbPoints();
    List<List<interpol_point>>& interpol_points = intpInfo.getIntPoints();

    getCurVelocity(interpol_points);
    List<label> intOrder = getIntOrder(interpol_points);

    forAll(interpol_points, ibp)
    {
        label cellI = cSurfCells[ibp];

        switch(intOrder[ibp])
        {
            case 0:
            {
                Ui[cellI] = ibPointsVal[ibp];
                break;
            }

            case 1:
            {
                vector VP1 = interpol_points[ibp][0].iVel_ - ibPointsVal[ibp];

                scalar deltaR = mag(interpol_points[ibp][0].iPoint_ - ibPoints[ibp]);
                scalar ds = mag(mesh.C()[cellI] - ibPoints[ibp]);

                vector linCoeff = VP1/(deltaR+SMALL);

                Ui[cellI] = linCoeff*ds + ibPointsVal[ibp];
                break;
            }

            case 2:
            {
                vector VP1 =  interpol_points[ibp][0].iVel_ - ibPointsVal[ibp];
                vector VP2 =  interpol_points[ibp][1].iVel_ - ibPointsVal[ibp];

                scalar deltaR1 = mag(interpol_points[ibp][0].iPoint_ - ibPoints[ibp]);
                scalar deltaR2 = mag(interpol_points[ibp][1].iPoint_
                        - interpol_points[ibp][0].iPoint_);

                scalar ds = mag(mesh.C()[cellI] - ibPoints[ibp]);

                vector quadCoeff = (VP2 - VP1)*deltaR1 - VP1*deltaR2;
                quadCoeff       /= (deltaR1*deltaR2*(deltaR1 + deltaR2)+SMALL);

                vector linCoeff  = (VP1-VP2)*Foam::pow(deltaR1,2.0);
                linCoeff        += 2.0*VP1*deltaR1*deltaR2;
                linCoeff        += VP1*Foam::pow(deltaR2,2.0);
                linCoeff        /= (deltaR1*deltaR2*(deltaR1 + deltaR2)+SMALL);

                Ui[cellI] = quadCoeff*ds*ds + linCoeff*ds + ibPointsVal[ibp];
                break;
            }
        }
    }
}
//---------------------------------------------------------------------------//
void line_interpol::getCurVelocity
(
    List<List<interpol_point>>& interpol_points
)
{
    List<DynamicPointList> interpol_pointToSync(Pstream::nProcs());
    List<DynamicLabelList> cellLabelToSync(Pstream::nProcs());
    List<DynamicList<Tuple2<label,label>>> indexesToSync(Pstream::nProcs());

    forAll(interpol_points, ibCellI)
    {
        forAll(interpol_points[ibCellI], iPoint)
        {
            interpol_point& curIPoint = interpol_points[ibCellI][iPoint];

            if(curIPoint.iProc_ == Pstream::myProcNo())
            {
                curIPoint.iVel_ =  interpV_->interpolate(
                    curIPoint.iPoint_,
                    curIPoint.iCell_
                );
            }
            else
            {
                if(curIPoint.iProc_ != -1)
                {
                    interpol_pointToSync[curIPoint.iProc_].append(curIPoint.iPoint_);
                    cellLabelToSync[curIPoint.iProc_].append(curIPoint.iCell_);
                    indexesToSync[curIPoint.iProc_].append(
                        Tuple2<label,label>(ibCellI, iPoint)
                    );
                }
            }
        }
    }

    PstreamBuffers pBufsIntP(Pstream::commsTypes::nonBlocking);
    PstreamBuffers pBufsIntC(Pstream::commsTypes::nonBlocking);

    for (label proci = 0; proci < Pstream::nProcs(); proci++)
    {
        if(proci != Pstream::myProcNo())
        {
            UOPstream sendIntP(proci, pBufsIntP);
            UOPstream sendIntC(proci, pBufsIntC);

            sendIntP << interpol_pointToSync[proci];
            sendIntC << cellLabelToSync[proci];
        }
    }

    pBufsIntP.finishedSends();
    pBufsIntC.finishedSends();

    List<DynamicPointList> intPRecv(Pstream::nProcs());
    List<DynamicLabelList> intCRecv(Pstream::nProcs());

    for (label proci = 0; proci < Pstream::nProcs(); proci++)
    {
        if (proci != Pstream::myProcNo())
        {
            UIPstream recvIntP(proci, pBufsIntP);
            UIPstream recvIntC(proci, pBufsIntC);

            DynamicPointList recIntP (recvIntP);
            DynamicLabelList recIntC (recvIntC);

            intPRecv[proci] = recIntP;
            intCRecv[proci] = recIntC;
        }
    }

    pBufsIntP.clear();
    pBufsIntC.clear();

    List<DynamicVectorList> intVelRtrn(Pstream::nProcs());

    forAll(intPRecv, proci)
    {
        forAll(intPRecv[proci], pi)
        {
            intVelRtrn[proci].append(
                interpV_->interpolate(
                    intPRecv[proci][pi],
                    intCRecv[proci][pi]
            ));
        }
    }

    for (label proci = 0; proci < Pstream::nProcs(); proci++)
    {
        if(proci != Pstream::myProcNo())
        {
            UOPstream sendIntV(proci, pBufsIntP);

            sendIntV << intVelRtrn[proci];
        }
    }

    pBufsIntP.finishedSends();

    List<DynamicVectorList> intVelRcv(Pstream::nProcs());

    for (label proci = 0; proci < Pstream::nProcs(); proci++)
    {
        if (proci != Pstream::myProcNo())
        {
            UIPstream recvIntV(proci, pBufsIntP);

            DynamicVectorList recIntV (recvIntV);

            intVelRcv[proci] = recIntV;
        }
    }

    pBufsIntP.clear();

    forAll(intVelRcv, proci)
    {
        forAll(intVelRcv[proci], pi)
        {
            Tuple2<label, label> cInd = indexesToSync[proci][pi];
            interpol_points[cInd.first()][cInd.second()].iVel_ = intVelRcv[proci][pi];
        }
    }
}
//---------------------------------------------------------------------------//
List<label> line_interpol::getIntOrder
(
    List<List<interpol_point>>& interpol_points
)
{
    List<label> intOrderToRtn(interpol_points.size(), 0);

    forAll(interpol_points, ibp)
    {
        forAll(interpol_points[ibp], ip)
        {
            if(interpol_points[ibp][ip].iProc_ != -1)
            {
                ++intOrderToRtn[ibp];
            }
        }
    }

    return intOrderToRtn;
}
//---------------------------------------------------------------------------//
