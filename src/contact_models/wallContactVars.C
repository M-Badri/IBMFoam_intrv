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
#include "wallContactVars.H"

using namespace Foam;

//---------------------------------------------------------------------------//
void wallContactVars::setMeanCntPars
(
    const fvMesh&   mesh,
    DynamicList<Tuple2<label,string>>& contactFaces,
    HashTable<physicalProperties,string,Hash<string>>& wallMeanPars
)
{
    scalar overallArea = 0;
    physicalProperties_.aY_ = 0;
    physicalProperties_.aG_ = 0;
    physicalProperties_.aMu_ = 0;
    physicalProperties_.maxAdhN_ = 0;
    physicalProperties_.reduceBeta_ = 0;

    forAll(contactFaces,faceI)
    {
        scalar area = mesh.magSf()[contactFaces[faceI].first()];
        overallArea += area;

        physicalProperties& cMeanCntPars(
            wallMeanPars[contactFaces[faceI].second()]
        );

        physicalProperties_.aY_ += (cMeanCntPars.aY_*area);
        physicalProperties_.aG_ += (cMeanCntPars.aG_*area);
        physicalProperties_.aMu_ += (cMeanCntPars.aMu_*area);
        physicalProperties_.maxAdhN_ += (cMeanCntPars.maxAdhN_*area);
        physicalProperties_.reduceBeta_ += (cMeanCntPars.reduceBeta_*area);
    }

    physicalProperties_.aY_ /= overallArea;
    physicalProperties_.aG_ /= overallArea;
    physicalProperties_.aMu_ /= overallArea;
    physicalProperties_.maxAdhN_ /= overallArea;
    physicalProperties_.reduceBeta_ /= overallArea;
}
//---------------------------------------------------------------------------//
void wallContactVars::setMeanCntPars_Plane
(
    List<scalar>& contactAreas,
    List<string> contactFaces,
    HashTable<physicalProperties,string,Hash<string>>& wallMeanPars
)
{
    scalar overallArea = 0;
    physicalProperties_.aY_ = 0;
    physicalProperties_.aG_ = 0;
    physicalProperties_.aMu_ = 0;
    physicalProperties_.maxAdhN_ = 0;
    physicalProperties_.reduceBeta_ = 0;

    forAll(contactFaces,faceI)
    {
        scalar area = contactAreas[faceI];
        overallArea += area;

        physicalProperties& cMeanCntPars(
            wallMeanPars[contactFaces[faceI]]
        );

        physicalProperties_.aY_ += (cMeanCntPars.aY_*area);
        physicalProperties_.aG_ += (cMeanCntPars.aG_*area);
        physicalProperties_.aMu_ += (cMeanCntPars.aMu_*area);
        physicalProperties_.maxAdhN_ += (cMeanCntPars.maxAdhN_*area);
        physicalProperties_.reduceBeta_ += (cMeanCntPars.reduceBeta_*area);
    }

    physicalProperties_.aY_ /= overallArea;
    physicalProperties_.aG_ /= overallArea;
    physicalProperties_.aMu_ /= overallArea;
    physicalProperties_.maxAdhN_ /= overallArea;
    physicalProperties_.reduceBeta_ /= overallArea;
}
//---------------------------------------------------------------------------//
