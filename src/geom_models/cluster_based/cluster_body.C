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
#include "cluster_body.H"

using namespace Foam;

//---------------------------------------------------------------------------//

void cluster_body::create_immersed_body
(
    volScalarField& body,
    Field<label>& octreeField,
    List<labelList>& cellPoints
)
{
    for(std::shared_ptr<geom_model>& gModel : ibGeomModelList)
    {
        gModel->create_immersed_body(
            body,
            octreeField,
            cellPoints
        );

        Info << "Periodic body created" << " mass: " << gModel->getM() << endl;
        Info << "Periodic body created" << " bbox: " << gModel->getBounds().min() << " " << gModel->getBounds().max() << endl;
        Info << "Periodic body created" << " intList: " << gModel->getInternalCellList()[Pstream::myProcNo()].size() << endl;
    }
}
//---------------------------------------------------------------------------//
void cluster_body::updateSurfList()
{
    surfCells_.clear();
    surfCells_.setSize(Pstream::nProcs());

    for(std::shared_ptr<geom_model>& gModel : ibGeomModelList)
    {
        surfCells_[Pstream::myProcNo()].append(
            gModel->getSurfaceCellList()[Pstream::myProcNo()]
        );
    }
}
//---------------------------------------------------------------------------//
void cluster_body::updateIntList()
{
    intCells_.clear();
    intCells_.setSize(Pstream::nProcs());

    for(std::shared_ptr<geom_model>& gModel : ibGeomModelList)
    {
        intCells_[Pstream::myProcNo()].append(
            gModel->getInternalCellList()[Pstream::myProcNo()]
        );
    }
}
//---------------------------------------------------------------------------//
void cluster_body::calculateGeometricalProperties
(
    volScalarField& body
)
{
    M_      = scalar(0);
    I_      = symmTensor::zero;

    for(std::shared_ptr<geom_model>& gModel : ibGeomModelList)
    {
        gModel->calculateGeometricalProperties(body);
        M_ += gModel->getM();
        I_ += gModel->getI();
    }
}
//---------------------------------------------------------------------------//
void cluster_body::calculateGeometricalPropertiesParallel
(
    volScalarField& body
)
{
    M_      = scalar(0);
    I_      = symmTensor::zero;

    for(std::shared_ptr<geom_model>& gModel : ibGeomModelList)
    {
        gModel->calculateGeometricalPropertiesParallel(body);

        Info << "Periodic body calculated" << " mass: " << gModel->getM() << endl;
    }
}
//---------------------------------------------------------------------------//
void cluster_body::setMassAndInertia()
{
    M_      = scalar(0);
    I_      = symmTensor::zero;

    for(std::shared_ptr<geom_model>& gModel : ibGeomModelList)
    {
        M_ += gModel->getM();
        I_ += gModel->getI();
    }
}
//---------------------------------------------------------------------------//
vector cluster_body::getCoM()
{
    return ibGeomModelList[0]->getCoM();
}
//---------------------------------------------------------------------------//
boundBox cluster_body::getBounds()
{
    DynamicPointList allBounds;
    for(std::shared_ptr<geom_model>& gModel : ibGeomModelList)
    {
        boundBox bBoxI = gModel->getBounds();
        allBounds.append(bBoxI.min());
        allBounds.append(bBoxI.max());
    }

    return boundBox(allBounds);
}
//---------------------------------------------------------------------------//
void cluster_body::synchronPos(label owner)
{
    for(std::shared_ptr<geom_model>& gModel : ibGeomModelList)
    {
        gModel->synchronPos(owner_);
    }
}
//---------------------------------------------------------------------------//
boolList cluster_body::pointInside(pointField pointF)
{
    boolList inside(pointF.size());

    forAll(pointF,pointI)
    {
        bool pointInside = false;
        for(std::shared_ptr<geom_model>& gModel : ibGeomModelList)
        {
            pointInside = gModel->pointInside(pointF[pointI]);
            if(pointInside)
                break;
        }
        inside[pointI] = pointInside;
    }

    return inside;
}
//---------------------------------------------------------------------------//
bool cluster_body::pointInside(point pointI)
{
    for(std::shared_ptr<geom_model>& gModel : ibGeomModelList)
    {
        if(gModel->pointInside(pointI))
            return true;
    }
    return false;
}
//---------------------------------------------------------------------------//
void cluster_body::getClosestPointAndNormal
(
    const point& startPoint,
    const vector& span,
    point& closestPoint,
    vector& normal
)
{
    List<point> closestPoints(ibGeomModelList.size());
    List<vector> closestNormals(ibGeomModelList.size());

    for(size_t ibI = 0; ibI < ibGeomModelList.size(); ++ibI)
    {
        ibGeomModelList[ibI]->getClosestPointAndNormal(
            startPoint,
            span,
            closestPoints[ibI],
            closestNormals[ibI]
        );
    }

    closestPoint = closestPoints[0];
    normal = closestNormals[0];
    for (int i = 1; i < closestPoints.size(); ++i)
    {
        if(mag(startPoint - closestPoint)
            > mag(startPoint - closestPoints[i]))
        {
            closestPoint = closestPoints[i];
            normal = closestNormals[i];
        }
    }
}
//---------------------------------------------------------------------------//
void cluster_body::resetBody(volScalarField& body)
{
    for(std::shared_ptr<geom_model>& gModel : ibGeomModelList)
    {
        gModel->resetBody(body);
    }
}
//---------------------------------------------------------------------------//
List<std::shared_ptr<boundBox>> cluster_body::getBBoxes()
{
    List<std::shared_ptr<boundBox>> retList;
    for(std::shared_ptr<geom_model>& gModel : ibGeomModelList)
    {
        List<std::shared_ptr<boundBox>> bBoxI = gModel->getBBoxes();
        retList.append(bBoxI);
    }
    return retList;
}
//---------------------------------------------------------------------------//
