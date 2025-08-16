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
    struct for interpolation point

SourceFiles
\*---------------------------------------------------------------------------*/

#ifndef interpol_point_H
#define interpol_point_H

namespace Foam
{
    struct interpol_point
    {
        point iPoint_ = vector::zero;
        label iCell_ = 0;
        label iProc_ = -1;
        vector iVel_ = vector::zero;

        interpol_point():
        iPoint_(vector::zero),
        iCell_(0),
        iProc_(-1)
        {};

        interpol_point(point iPoint, label iCell, label iProc):
        iPoint_(iPoint),
        iCell_(iCell),
        iProc_(iProc)
        {};

        ~interpol_point(){};
    };
}
#endif
