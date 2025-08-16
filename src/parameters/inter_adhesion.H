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
    Static class to store inter-particle adhesion coefficients
\*---------------------------------------------------------------------------*/

#ifndef inter_adhesion_H
#define inter_adhesion_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{


class inter_adhesion
{
    public:

    static const HashTable<scalar,string,Hash<string>>& getInterAdhesion()
    {
        return matInterAdh_;
    };

    static void inter_adhesion_insert
    (
        const string& key,
        const scalar& value
    )
    {
        matInterAdh_.insert(key,value);
    };

    static HashTable<scalar,string,Hash<string>> matInterAdh_;

};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
