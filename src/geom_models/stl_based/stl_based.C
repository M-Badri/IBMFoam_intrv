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
#include "stl_based.H"

using namespace Foam;

//---------------------------------------------------------------------------//
stl_based::stl_based
(
    const  fvMesh&   mesh,
    const contactType cType,
    word      stlPath,
    scalar  thrSurf
)
:
geom_model(mesh,cType,thrSurf),
bodySurfMesh_
(
    IOobject
    (
        stlPath,
        mesh,
        IOobject::MUST_READ,
        IOobject::NO_WRITE
    )
),
stlPath_(stlPath)
{
    historyPoints_ = bodySurfMesh_.points();
    triSurf_.set(new triSurface(bodySurfMesh_));
    triSurfSearch_.set(new triSurfaceSearch(triSurf_()));
}
//---------------------------------------------------------------------------//
vector stl_based::add_modelReturnRandomPosition
(
    const bool allActiveCellsInMesh,
    const boundBox  cellZoneBounds,
    Random&          randGen
)
{
    vector ranVec(vector::zero);

    pointField bSMeshPts = bodySurfMesh_.points();

    vector CoM(vector::zero);
    forAll(bSMeshPts,point)
    {
        CoM += bSMeshPts[point];
    }
    CoM/= bSMeshPts.size();

    const vector validDirs = (geometricD + vector::one)/2;
    vector dirCorr(cmptMultiply((vector::one - validDirs),CoM));
    dirCorr += cmptMultiply((vector::one - validDirs),0.5*(mesh_.bounds().max() + mesh_.bounds().min()));

    boundBox bodySurfBounds(bSMeshPts);

    vector maxScales(cellZoneBounds.max() - bodySurfBounds.max());
    maxScales -= cellZoneBounds.min() - bodySurfBounds.min();
    maxScales *= 0.5*0.9;

    InfoH << add_model_Info << "-- add_modelMessage-- "
        << "acceptable movements: " << maxScales << endl;

    scalar ranNum = 0;
    for (int i=0;i<3;i++)
    {
        ranNum = 2.0*maxScales[i]*randGen.scalar01() - 1.0*maxScales[i];
        ranVec[i] = ranNum;
    }

    ranVec = cmptMultiply(validDirs,ranVec);                            
    ranVec += dirCorr;

    return ranVec;
}
//---------------------------------------------------------------------------//
void stl_based::bodyMovePoints
(
    vector translVec
)
{
    pointField bodyPoints(bodySurfMesh_.points());
    bodyPoints += translVec;

    bodySurfMesh_.movePoints(bodyPoints);
    triSurf_.reset(new triSurface(bodySurfMesh_));
    triSurfSearch_.reset(new triSurfaceSearch(triSurf_()));
}
//---------------------------------------------------------------------------//
void stl_based::bodyScalePoints
(
    scalar scaleFac
)
{
    pointField bodyPoints(bodySurfMesh_.points());

    vector CoM(vector::zero);
    forAll(bodyPoints,point)
    {
        CoM += bodyPoints[point];
    }
    CoM/= bodyPoints.size();

    bodyPoints -= CoM;
    bodySurfMesh_.movePoints(bodyPoints);
    bodySurfMesh_.scalePoints(scaleFac);
    bodyPoints = bodySurfMesh_.points();
    bodyPoints += CoM;
    bodySurfMesh_.movePoints(bodyPoints);
    triSurf_.reset(new triSurface(bodySurfMesh_));
    triSurfSearch_.reset(new triSurfaceSearch(triSurf_()));
}
//---------------------------------------------------------------------------//
void stl_based::bodyRotatePoints
(
    scalar rotAngle,
    vector axisOfRot
)
{
    pointField bodyPoints(bodySurfMesh_.points());
    vector CoM(vector::zero);
    forAll(bodyPoints,point)
    {
        CoM += bodyPoints[point];
    }
    CoM/= bodyPoints.size();

    tensor rotMatrix(Foam::cos(rotAngle)*tensor::I);

    rotMatrix += Foam::sin(rotAngle)*tensor(
            0.0,      -axisOfRot.z(),  axisOfRot.y(),
            axisOfRot.z(), 0.0,       -axisOfRot.x(),
        -axisOfRot.y(), axisOfRot.x(),  0.0
    );

    rotMatrix += (1.0-Foam::cos(rotAngle))*(axisOfRot * axisOfRot);

    bodyPoints -= CoM;
    bodyPoints = rotMatrix & bodyPoints;
    bodyPoints += CoM;
    bodySurfMesh_.movePoints(bodyPoints);
    triSurf_.reset(new triSurface(bodySurfMesh_));
    triSurfSearch_.reset(new triSurfaceSearch(triSurf_()));
}
//---------------------------------------------------------------------------//
void stl_based::synchronPos(label owner)
{
    PstreamBuffers pBufs(Pstream::commsTypes::nonBlocking);

    owner = (owner == -1) ? owner_ : owner;

    if (owner == Pstream::myProcNo())
    {
        for (label proci = 0; proci < Pstream::nProcs(); proci++)
        {
            UOPstream send(proci, pBufs);
            send << bodySurfMesh_.points();
        }
    }

    pBufs.finishedSends();
    UIPstream recv(owner, pBufs);
    pointField bodyPoints (recv);

    // move mesh
    bodySurfMesh_.movePoints(bodyPoints);
    triSurf_.reset(new triSurface(bodySurfMesh_));
    triSurfSearch_.reset(new triSurfaceSearch(triSurf_()));
}
//---------------------------------------------------------------------------//
volumeType stl_based::getVolumeType(sub_volume& sv, bool cIb)
{
    const indexedOctree<treeDataTriSurface>& tree = triSurfSearch_->tree();
    const treeDataTriSurface& shapes = tree.shapes();

    autoPtr<labelList>& shapesIn = sv.getVolumeInfo(cIb).shapesIn_;

    if (shapesIn.empty())
    {
        std::shared_ptr<sub_volume> parentSV = sv.parentSV();
        if (parentSV)
        {
            const labelList& parentShapesIn = *(parentSV->getVolumeInfo(cIb).shapesIn_);
            labelHashSet shapesInSV;

            forAll(parentShapesIn, i)
            {
                label shapeI = parentShapesIn[i];
                if (shapes.overlaps(shapeI, sv))
                {
                    shapesInSV.insert(shapeI);
                }
            }

            shapesIn.reset(new labelList(shapesInSV.toc()));
        }
        else
        {
            shapesIn.reset(new labelList(tree.findBox(sv)));
        }
    }

    if (shapesIn->size() > 0)
    {
        return volumeType::mixed;
    }

    return pointInside(sv.midpoint())
            ? volumeType::inside : volumeType::outside;
}
//---------------------------------------------------------------------------//
bool stl_based::limitFinalSubVolume
(
    const sub_volume& sv,
    bool cIb,
    boundBox& limBBox
)
{
    limBBox = boundBox(sv.min(), sv.max());
    const autoPtr<labelList>& shapesIn = sv.getVolumeInfo(cIb).shapesIn_;
    if (shapesIn.empty())
    {
        return false;
    }

    vector normal = vector::zero;
    DynamicPointList intersectionPoints;
    forAll(*shapesIn, i)
    {
        getIntersectionPoints((*shapesIn)[i], sv, intersectionPoints);
        normal += (*triSurf_)[(*shapesIn)[i]].area(triSurf_->points());
    }

    if (intersectionPoints.size() == 0)
    {
        return false;
    }

    point closestPoint = vector::zero;
    forAll(intersectionPoints, i)
    {
        closestPoint += intersectionPoints[i];
    }
    closestPoint /= intersectionPoints.size();

    scalar magX = mag(normal.x());
    scalar magY = mag(normal.y());
    scalar magZ = mag(normal.z());

    if (magX > magY && magX > magZ)
    {
        if (sign(normal.x()) < 0)
        {
            limBBox.min().x() = closestPoint.x();
            return true;
        }
        else
        {
            limBBox.max().x() = closestPoint.x();
            return true;
        }
    }
    else if (magY > magX && magY > magZ)
    {
        if (sign(normal.y()) < 0)
        {
            limBBox.min().y() = closestPoint.y();
            return true;
        }
        else
        {
            limBBox.max().y() = closestPoint.y();
            return true;
        }
    }
    else if (magZ > magX && magZ > magY)
    {
        if (sign(normal.z()) < 0)
        {
            limBBox.min().z() = closestPoint.z();
            return true;
        }
        else
        {
            limBBox.max().z() = closestPoint.z();
            return true;
        }
    }

    return false;
    
}
//---------------------------------------------------------------------------//
void stl_based::getIntersectionPoints
(
    const label index,
    const treeBoundBox& cubeBb,
    DynamicPointList& intersectionPoints
)
{
    const pointField& points = triSurf_->points();
    const typename triSurface::FaceType& f = (*triSurf_)[index];

    for (auto ind : f)
    {
        if (cubeBb.contains(points[ind]))
        {
            intersectionPoints.append(points[ind]);
        }
    }

    const point fc = f.centre(points);

    if (f.size() == 3)
    {
        return intersectBb
        (
            points[f[0]],
            points[f[1]],
            points[f[2]],
            cubeBb,
            intersectionPoints
        );
    }
    else
    {
        forAll(f, fp)
        {
            intersectBb
            (
                points[f[fp]],
                points[f[f.fcIndex(fp)]],
                fc,
                cubeBb,
                intersectionPoints
            );
        }
    }

    return;
}
//---------------------------------------------------------------------------//
void stl_based::intersectBb
(
    const point& p0,
    const point& p1,
    const point& p2,
    const treeBoundBox& cubeBb,
    DynamicPointList& intersectionPoints
)
{
    const vector p10 = p1 - p0;
    const vector p20 = p2 - p0;

    const point& min = cubeBb.min();
    const point& max = cubeBb.max();

    const point& cube0 = min;
    const point  cube1(min.x(), min.y(), max.z());
    const point  cube2(max.x(), min.y(), max.z());
    const point  cube3(max.x(), min.y(), min.z());

    const point  cube4(min.x(), max.y(), min.z());
    const point  cube5(min.x(), max.y(), max.z());
    const point  cube7(max.x(), max.y(), min.z());


    point pInter;
    pointField origin(4);

    origin[0] = cube0;
    origin[1] = cube1;
    origin[2] = cube5;
    origin[3] = cube4;

    scalar maxSx = max.x() - min.x();

    if (triangleFuncs::intersectAxesBundle(p0, p10, p20, 0, origin, maxSx, pInter))
    {
        intersectionPoints.append(pInter);
    }


    origin[0] = cube0;
    origin[1] = cube1;
    origin[2] = cube2;
    origin[3] = cube3;

    scalar maxSy = max.y() - min.y();

    if (triangleFuncs::intersectAxesBundle(p0, p10, p20, 1, origin, maxSy, pInter))
    {
        intersectionPoints.append(pInter);
    }


    origin[0] = cube0;
    origin[1] = cube3;
    origin[2] = cube7;
    origin[3] = cube4;

    scalar maxSz = max.z() - min.z();

    if (triangleFuncs::intersectAxesBundle(p0, p10, p20, 2, origin, maxSz, pInter))
    {
        intersectionPoints.append(pInter);
    }



    if (cubeBb.intersects(p0, p1, pInter))
    {
        intersectionPoints.append(pInter);
    }
    if (cubeBb.intersects(p1, p2, pInter))
    {
        intersectionPoints.append(pInter);
    }
    if (cubeBb.intersects(p2, p0, pInter))
    {
        intersectionPoints.append(pInter);
    }
}
//---------------------------------------------------------------------------//
void stl_based::setBodyPosition(pointField pos)
{
    bodySurfMesh_.movePoints(pos);
    triSurf_.reset(new triSurface(bodySurfMesh_));
    triSurfSearch_.reset(new triSurfaceSearch(triSurf_()));
}
//---------------------------------------------------------------------------//
