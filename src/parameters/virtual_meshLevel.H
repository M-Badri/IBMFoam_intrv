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

#ifndef virtual_meshLevel_H
#define virtual_meshLevel_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{



class virtual_meshLevel
{
    public:

    static void setVirtualMeshLevel(scalar virtual_meshLevel,scalar charCellSize)
    {
        virtual_meshLevel_ = virtual_meshLevel;
        levelOfDivision_ = int(pow(2,virtual_meshLevel_));
        charCellSize_ = charCellSize;
    }

    static scalar getVirtualMeshLevel()
    {
        return virtual_meshLevel_;
    }

    static scalar getLevelOfDivision()
    {
        return levelOfDivision_;
    }

    static scalar getCharCellSize()
    {
        return charCellSize_;
    }

    static scalar virtual_meshLevel_;
    static scalar levelOfDivision_;
    static scalar charCellSize_;

};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
