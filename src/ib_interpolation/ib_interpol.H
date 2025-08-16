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
    ib_interpol.C
\*---------------------------------------------------------------------------*/

#ifndef ib_interpol_H
#define ib_interpol_H

#include "interpol_info.H"

namespace Foam
{
    class ib_interpol
    {
        typedef DynamicList<label> DynamicLabelList;
        typedef DynamicList<scalar> DynamicScalarList;
        typedef DynamicList<vector> DynamicVectorList;
        typedef DynamicList<point> DynamicPointList;

        public:

        // Constructors ------------------------------------------------
        ib_interpol(){};

        // Destructors -------------------------------------------------
        virtual ~ib_interpol(){};

        virtual void resetInterpolator(volVectorField& V){};

        virtual void ib_interpolate(
            interpol_info& intpInfo,
            volVectorField& Ui,
            vectorField ibPointsVal,
            const Foam::fvMesh& mesh
        ) = 0;

        virtual interpol_info* createIntpInfo
        (
            const Foam::fvMesh& mesh,
            std::shared_ptr<geom_model>& gModel
        ) = 0;
    };
}
#endif
