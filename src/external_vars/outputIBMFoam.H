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
    methods to solve particle-particle contact

SourceFiles
    prtContact.C
\*---------------------------------------------------------------------------*/

#ifndef outputIBMFoam_H
#define outputIBMFoam_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{



namespace outputIBMFoam
{
    enum outputType
    {
        basic_Info,
        iB_Info,
        DEM_Info,
        add_model_Info,
        parallelDEM_Info
    };

    class
    InfoIBMFoam
    {
        private:
        bool printBasic;
        bool printIb;
        bool printDEM;
        bool printAddModel;
        bool printParallelDEM;
        bool print;
        bool printParDEM;

        HashTable<bool, outputType,Hash<label>> outputTable;

        public:
        // Constructors ----------------------------------------------------
        InfoIBMFoam()
        {
            printBasic = true;
            printIb = true;
            printDEM = true;
            printAddModel = true;
            printParallelDEM = true;
            print = true;
            printParDEM = true;

            outputTable.insert(basic_Info, printBasic);
            outputTable.insert(iB_Info, printIb);
            outputTable.insert(DEM_Info, printDEM);
            outputTable.insert(add_model_Info, printAddModel);
            outputTable.insert(parallelDEM_Info, printParallelDEM);
        }
        // Destructors -----------------------------------------------------
        ~InfoIBMFoam(){}

        void setOutput(bool basicB, bool ibB, bool demB, bool add_modelB, bool parDEMB)
        {
            outputTable[basic_Info] = basicB;
            outputTable[iB_Info] = ibB;
            outputTable[DEM_Info] = demB;
            outputTable[add_model_Info] = add_modelB;
            outputTable[parallelDEM_Info] = parDEMB;
        }

        InfoIBMFoam& operator<<(const outputType& data)
        {
            if(outputTable[data]  && data != parallelDEM_Info )
            {
                print = true;
                printParDEM = false;
            }
            else if (outputTable[data] && data == parallelDEM_Info)
            {
                printParDEM = true;
                print = false;
            }
            else
            {
                print = false;
                printParDEM = false;
            }
            return *this;
        }

        template<typename T>
        InfoIBMFoam& operator<<(const T& data)
        {
            if(printParDEM && !print)
            {
                Pout << data;
            }
            if(print && !printParDEM)
            {
                Info << data;
            }
            return *this;
        }

        typedef ostream& (manip)(ostream&);
        InfoIBMFoam& operator<<(manip m)
        {
            if(print && !printParDEM)
            {
                Info << m;
            }

            if(printParDEM && !print)
            {
                Pout << m;
            }
            return *this;
        }
    };
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

#include "dec_externVars.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
