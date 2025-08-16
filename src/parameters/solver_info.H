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
    contactModel

Description
    Static class for storing virtual mesh level
\*---------------------------------------------------------------------------*/

#ifndef solver_info_H
#define solver_info_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{



class solver_info
{
    public:

    static void setOnlyDEM(bool onlyDEM)
    {
        onlyDEM_ = onlyDEM;
    }

    static void setNSolidsTreshnold(label nSolidsTreshold)
    {
        nSolidsTreshnold_ = nSolidsTreshold;
    }

    static bool getOnlyDEM()
    {
        return onlyDEM_;
    }

    static label getNSolidsTreshnold()
    {
        return nSolidsTreshnold_;
    }

    static bool onlyDEM_;
    static label nSolidsTreshnold_;

};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
