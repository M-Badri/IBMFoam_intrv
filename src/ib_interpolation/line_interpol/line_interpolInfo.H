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

#ifndef line_interpolInfo_H
#define line_interpolInfo_H

#include "interpol_info.H"
#include "interpol_point.H"

namespace Foam
{
    class line_interpolInfo
    : public interpol_info
    {
        typedef DynamicList<label> DynamicLabelList;
        typedef DynamicList<scalar> DynamicScalarList;
        typedef DynamicList<vector> DynamicVectorList;
        typedef DynamicList<point> DynamicPointList;

        private:

        autoPtr<List<List<interpol_point>>> interpol_pointsPtr_; // Points used for interp.
        autoPtr<List<List<vector>>> intVelocPtr_; // Velocity at int points

        interpol_point  findIntPoint
        (
            interpol_point& fromP,
            point& endP
        );

        void correctIntPoint
        (
            point ibPoint,
            interpol_point& cPoint
        );

        vector getClosestPoint
        (
            vector ibPoint,
            interpol_point& cPoint
        );

        label getFaceInDir
        (
            const interpol_point& retPoint,
            const label prevFace
        );

        bool pointInCell
        (
            point pToCheck,
            label cToCheck
        );

        void syncIntPoints();

        public:

        // Constructors ------------------------------------------------
        line_interpolInfo
        (
            const  fvMesh&   mesh,
            std::shared_ptr<geom_model>& gModel
        );

        // Destructors -------------------------------------------------
        ~line_interpolInfo();

        void setIntpInfo() final override;

        void resetIntpInfo(label newSize)
        {
            ibPointsPtr_.reset(new List<point>(newSize, vector::zero));
            ibNormalsPtr_.reset(new List<point>(newSize, vector::zero));
            interpol_pointsPtr_.reset(new List<List<interpol_point>>(newSize));
            intVelocPtr_.reset(new List<List<vector>>(newSize));
        }

        List<List<interpol_point>>& getIntPoints()
        {
            return interpol_pointsPtr_();
        }

        const List<List<interpol_point>>& getIntPoints() const
        {
            return interpol_pointsPtr_();
        }

        List<List<vector>>& getIntVeloc()
        {
            return intVelocPtr_();
        }

        const List<List<vector>>& getIntVeloc() const
        {
            return intVelocPtr_();
        }
    };
}
#endif
