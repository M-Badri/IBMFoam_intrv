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
    DEM parameters
\*---------------------------------------------------------------------------*/

#ifndef parameters_H
#define parameters_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "materialProperties.H"
#include "inter_adhesion.H"
#include "wall_mat_info.H"
#include "wall_plane_info.H"
#include "virtual_meshLevel.H"
#include "cyclic_plane_info.H"
#include "contact_model_info.H"
#include "solver_info.H"

namespace Foam
{


HashTable<material_info,string,Hash<string>> materialProperties::materialProps_ =
    HashTable<material_info,string,Hash<string>>();

HashTable<scalar,string,Hash<string>> inter_adhesion::matInterAdh_ =
    HashTable<scalar,string,Hash<string>>();

HashTable<material_info,string,Hash<string>> wall_mat_info::wall_mat_info_ =
    HashTable<material_info,string,Hash<string>>();

HashTable<List<vector>,string,Hash<string>> wall_plane_info::wall_plane_info_ =
    HashTable<List<vector>,string,Hash<string>>();

HashTable<List<vector>,string,Hash<string>> cyclic_plane_info::cyclic_plane_info_ =
    HashTable<List<vector>,string,Hash<string>>();

std::map<string, string> cyclic_plane_info::cyclicPatchPairs_ =
    std::map<string, string>();

scalar virtual_meshLevel::virtual_meshLevel_ = 1;
scalar virtual_meshLevel::levelOfDivision_ = 1;
scalar virtual_meshLevel::charCellSize_ = 1;
scalar contact_model_info::LcCoeff_ = 4.0;
label  solver_info::nSolidsTreshnold_ = 1000;
bool solver_info::onlyDEM_ = false;
bool contact_model_info::useChenRotationalModel_ = false;
bool contact_model_info::useMindlinRotationalModel_ = false;
bool contact_model_info::increasedDamping_ = false;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //