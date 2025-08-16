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
#include "LS_interpol.H"
#include "scalarMatrices.H"

using namespace Foam;

//---------------------------------------------------------------------------//
LS_interpol::LS_interpol
(
    const fvMesh& mesh,
    scalar distFactor,
    scalar radiusFactor,
    scalar angleFactor,
    scalar maxCCRows
)
:
mesh_(mesh),
distFactor_(distFactor),
radiusFactor_(radiusFactor),
angleFactor_(angleFactor),
maxCCRows_(maxCCRows)
{}
LS_interpol::~LS_interpol()
{}
//---------------------------------------------------------------------------//
void LS_interpol::ib_interpolate
(
    interpol_info& intpInfo,
    volVectorField& Ui,
    vectorField ibPointsVal,
    const Foam::fvMesh& mesh
)
{
    LS_interpolInfo& lsInfo
        = dynamic_cast<LS_interpolInfo&>(intpInfo);

    imposeDirichletCondition
    (
        lsInfo,
        Ui,
        ibPointsVal,
        mesh
    );
}
//---------------------------------------------------------------------------//
void LS_interpol::imposeDirichletCondition
(
    LS_interpolInfo& intpInfo,
    volVectorField& Ui,
    vectorField& ibPointsVal,
    const Foam::fvMesh& mesh
)
{
    const DynamicLabelList& ibCells = intpInfo.getSurfCells();
    const List<point>& ibPoints = intpInfo.getIbPoints();
    const List<vector>& ibNormals = intpInfo.getIbNormals();
    const labelListList& cellCells = intpInfo.getCellCells();
    const PtrList<scalarRectangularMatrix>& cInvDirMats =
        intpInfo.getInvDirMats();

    autoPtr<vectorField> tpolyU(new vectorField(Ui.internalField(),ibCells));
    vectorField& polyU = tpolyU();

    const vectorField& C = mesh_.cellCentres();

    label nCoeffs = 5;
    if(case3D)
    {
        nCoeffs += 4;
    }

    label counter = 0;
    scalarField error(ibCells.size(), 0);
    scalar maxError = 0;

    do
    {
        counter++;
        scalar maxMagPolyU = 0;
        forAll(ibCells, cellI)
        {
            label curCell = ibCells[cellI];
            const labelList& curCells = cellCells[cellI];
            const scalarRectangularMatrix& curInvMatrix = cInvDirMats[cellI];

            if(curCells.size() < nCoeffs)
            {
                Info << "Not eneough interp Points" << endl;
            }
            else
            {
                vectorField coeffs(nCoeffs, vector::zero);
                vectorField source(curCells.size(), vector::zero);

                for(label i = 0; i < curCells.size(); i++)
                {
                    source[i] = Ui[curCells[i]] - ibPointsVal[cellI];
                }

                for(label i = 0; i < nCoeffs; i++)
                {
                    for(label j = 0; j < source.size(); j++)
                    {
                        coeffs[i] += curInvMatrix[i][j]*source[j];
                    }
                }

                vector oldPolyU = polyU[cellI];
                vector R = C[curCell] - ibPoints[cellI];
                if((R & ibNormals[cellI]) < 0)
                    R = vector::zero;

                if(case3D)
                {
                    polyU[cellI] = ibPointsVal[cellI]
                        + coeffs[0]*R.x()
                        + coeffs[1]*R.y()
                        + coeffs[2]*R.x()*R.y()
                        + coeffs[3]*sqr(R.x())
                        + coeffs[4]*sqr(R.y())
                        + coeffs[5]*R.z()
                        + coeffs[6]*R.x()*R.z()
                        + coeffs[7]*R.y()*R.z()
                        + coeffs[8]*sqr(R.z());
                }
                else
                {
                    labelList A (2,0.0);
                    label validDim = 0;
                    for(label dim = 0; dim < 3; dim++)
                    {
                        if(dim != emptyDim)
                        {
                            A[validDim++] = dim;
                        }
                    }
                    polyU[cellI] = ibPointsVal[cellI]
                        + coeffs[0]*R[A[0]]
                        + coeffs[1]*R[A[1]]
                        + coeffs[2]*R[A[0]]*R[A[1]]
                        + coeffs[3]*sqr(R[A[0]])
                        + coeffs[4]*sqr(R[A[1]]);
                }

                error[cellI] = mag(polyU[cellI] - oldPolyU);
            }
        }

        forAll(ibCells, ci)
        {
            Ui[ibCells[ci]] = polyU[ci];
        }

        if(!polyU.empty())
        {
            maxMagPolyU = max(mag(polyU));
        }
        else
        {
            maxMagPolyU = 0;
        }

        error /= maxMagPolyU + SMALL;

        if(!polyU.empty())
        {
            maxError = max(error);
        }
        else
        {
            maxError = 0;
        }

    } while(maxError > 0.01 && counter < 5);
}
//---------------------------------------------------------------------------//
void LS_interpol::adjustPhi
(
    LS_interpolInfo& intpInfo,
    surfaceScalarField& phi
)
{

    scalar fluxIn = 0;
    scalar fluxOut = 0;

    const labelList& ibFaces = intpInfo.getIbFaces();
    const boolList& ibFaceFlips = intpInfo.getIbFacesFlips();

    forAll (ibFaces, faceI)
    {
        const label curFace = ibFaces[faceI];
        const bool curFlip = ibFaceFlips[faceI];

        const scalar curFlux = phi[curFace];

        if (!curFlip)
        {
            if (curFlux >= 0)
            {
                fluxOut += curFlux;
            }
            else
            {
                fluxIn -= curFlux;
            }
        }
        else
        {
            if (curFlux >= 0)
            {
                fluxIn += curFlux;
            }
            else
            {
                fluxOut -= curFlux;
            }
        }
    }

    scalar imbalance = fluxIn - fluxOut;

    Info<< " fluxIn = " << fluxIn << " fluxOut = " << fluxOut
        << " imbalance = " << imbalance
        << endl;

    scalar massCorr = 1.0;

    if (mag(imbalance) > SMALL)
    {
        if (fluxIn > fluxOut)
        {
            massCorr = 1 - imbalance/(fluxIn + SMALL);

            Info<< "Scaling down incoming flux with factor = "
                << massCorr << endl;

            scalar newFluxIn = 0;

            forAll (ibFaces, faceI)
            {
                const label curFace = ibFaces[faceI];
                const bool curFlip = ibFaceFlips[faceI];

                if (mesh_.isInternalFace(curFace))
                {
                    scalar& curFlux = phi[curFace];

                    if (!curFlip)
                    {
                        if (curFlux < 0)
                        {
                            curFlux *= massCorr;
                            newFluxIn -= curFlux;
                        }
                    }
                    else
                    {
                        if (curFlux >= 0)
                        {
                            curFlux *= massCorr;
                            newFluxIn += curFlux;
                        }
                    }
                }
            }
        }
        else
        {
            massCorr = 1 + imbalance/(fluxOut + SMALL);

            Info<< "Scaling down outgoing flux with factor = "
                << massCorr << endl;

            scalar newFluxOut = 0;

            forAll (ibFaces, faceI)
            {
                const label curFace = ibFaces[faceI];
                const bool curFlip = ibFaceFlips[faceI];

                if (mesh_.isInternalFace(curFace))
                {
                    scalar& curFlux = phi[curFace];

                    if (!curFlip)
                    {
                        if (curFlux >= 0)
                        {
                            curFlux *= massCorr;
                            newFluxOut += curFlux;
                        }
                    }
                    else
                    {
                        if (curFlux < 0)
                        {
                            curFlux *= massCorr;
                            newFluxOut -= curFlux;
                        }
                    }
                }
            }
        }
    }
}
//---------------------------------------------------------------------------//
