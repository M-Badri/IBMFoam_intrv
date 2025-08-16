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

#ifndef contact_model_info_H
#define contact_model_info_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{


class contact_model_info
{
    public:

    static void setLcCoeff(scalar LcCoeff)
    {
        LcCoeff_ = LcCoeff;
    }

    static void set_rotation_model(label modelID)
    {
        if (modelID == 0)
        {
            useChenRotationalModel_ = true;
            useMindlinRotationalModel_ = false;
        }
        else if (modelID == 1)
        {
            useChenRotationalModel_ = false;
            useMindlinRotationalModel_ = true;
        }
    }

    static void set_increased_damping(bool increasedDamping)
    {
        increasedDamping_ = increasedDamping;
    }

    static scalar getLcCoeff()
    {
        return LcCoeff_;
    }

    static bool getUseChenRotationalModel()
    {
        return useChenRotationalModel_;
    }

    static bool getUseMindlinRotationalModel()
    {
        return useMindlinRotationalModel_;
    }

    static bool getIncreasedDamping()
    {
        return increasedDamping_;
    }

    static scalar LcCoeff_;
    static bool useChenRotationalModel_;
    static bool useMindlinRotationalModel_;
    static bool increasedDamping_;



};


} // End namespace Foam


#endif

// ************************************************************************* //
