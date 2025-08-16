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

Description
    class for interpolation

SourceFiles
\*---------------------------------------------------------------------------*/

#ifndef LS_interpolInfo_H
#define LS_interpolInfo_H

#include "interpol_info.H"

namespace Foam
{
    class LS_interpolInfo
    : public interpol_info
    {
        typedef DynamicList<label> DynamicLabelList;
        typedef DynamicList<scalar> DynamicScalarList;
        typedef DynamicList<vector> DynamicVectorList;
        typedef DynamicList<point> DynamicPointList;

        private:

        autoPtr<labelListList> cellCellsPtr_; 
        PtrList<scalarRectangularMatrix> invDirichletMatrices_;
        autoPtr<DynamicLabelList> ibFacesPtr_;
        autoPtr<DynamicList<bool>> ibFaceFlipsPtr_;

        const  scalar   distFactor_;
        const  scalar   radiusFactor_;
        const  scalar   angleFactor_;
        const  scalar   maxCCRows_;

        void findCellCells
        (
            const label cellId,
            labelList& cellCells,
            scalar& centerMeanDist
        );

        void getInvDirichletMatrix
        (
            const DynamicLabelList& ibCells
        );

        void getIbCellsFaces
        (
            const DynamicLabelList& ibCells
        );

        scalarSquareMatrix LUinvert
        (
            scalarSquareMatrix& matrix
        );


        public:

        // Constructors ------------------------------------------------
        LS_interpolInfo
        (
            const  fvMesh&   mesh,
            std::shared_ptr<geom_model>& gModel,
            scalar distFactor,
            scalar radiusFactor,
            scalar angleFactor,
            scalar maxCCRows
        );

        // Destructors -------------------------------------------------
        ~LS_interpolInfo();

        void resetIntpInfo(label newSize)
        {
            ibPointsPtr_.reset(new List<point>(newSize, vector::zero));
            ibNormalsPtr_.reset(new List<vector>(newSize, vector::zero));
            cellCellsPtr_.reset(new labelListList(newSize));
            ibFacesPtr_.reset(new DynamicLabelList());
            ibFaceFlipsPtr_.reset(new DynamicList<bool>());
        }

        void resetInvDirMat(label newSize)
        {
            invDirichletMatrices_.setSize(newSize);
        }

        void setIntpInfo() final override;

        List<vector>& getIbNormals()
        {
            return ibNormalsPtr_();
        }

        const List<vector>& getIbNormals() const
        {
            return ibNormalsPtr_();
        }

        labelListList& getCellCells()
        {
            return cellCellsPtr_();
        }

        const labelListList& getCellCells() const
        {
            return cellCellsPtr_();
        }

        PtrList<scalarRectangularMatrix>& getInvDirMats()
        {
            return invDirichletMatrices_;
        }

        const PtrList<scalarRectangularMatrix>& getInvDirMats() const
        {
            return invDirichletMatrices_;
        }

        DynamicLabelList& getIbFaces()
        {
            return ibFacesPtr_();
        }

        const DynamicLabelList& getIbFaces() const
        {
            return ibFacesPtr_();
        }

        DynamicList<bool>& getIbFacesFlips()
        {
            return ibFaceFlipsPtr_();
        }

        const DynamicList<bool>& getIbFacesFlips() const
        {
            return ibFaceFlipsPtr_();
        }
    };
}
#endif
