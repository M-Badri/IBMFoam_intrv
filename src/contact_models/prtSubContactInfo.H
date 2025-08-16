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
    class for storage of contact variables

SourceFiles
    prtSubContactInfo.C
\*---------------------------------------------------------------------------*/

#ifndef prtSubContactInfo_H
#define prtSubContactInfo_H

#include "dictionary.H"
#include "fvCFD.H"

#include "geom_model.H"
#include "ibContactClass.H"
#include "ibContactVars.H"
#include "prtContactVars.H"
#include "virtual_mesh.H"
#include "physicalProperties.H"
#include "forces.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{


class prtSubContactInfo
{
    private:

    typedef DynamicList<label> DynamicLabelList;
    typedef DynamicList<scalar> DynamicScalarList;
    typedef DynamicList<vector> DynamicVectorList;
    typedef DynamicList<point> DynamicPointList;

    Tuple2<label,label> contactPair_;

    physicalProperties physicalProperties_;

    std::shared_ptr<virtual_meshInfo> vmInfo_;

    prtContactVars prtCntVars_;

    vector cLVec_;       
    vector tLVec_;       

    vector cVeli_;       
    vector tVeli_;       

    scalar Vn_;          
    scalar Lc_;          

    Tuple2<forces,forces> outForce_; 

    vector FtPrev_ = vector::zero;     

    bool contact_resolved_ = false;

    vector getVeli(ibContactVars& cVars, vector& lVec);

    public:

    // Constructors ------------------------------------------------
    prtSubContactInfo
    (
        const Tuple2<label,label>& contactPair,
        const physicalProperties& physicalProperties
    );

    // Destructors -------------------------------------------------
    ~prtSubContactInfo();

    const Tuple2<label, label>& getCPair()
    {
        return contactPair_;
    }

    physicalProperties& getPhysicalProperties()
    {
        return physicalProperties_;
    }

    prtContactVars& getprtCntVars()
    {
        return prtCntVars_;
    }

    void evalVariables(
        ibContactClass& cCoM,
        ibContactClass& tCoM,
        ibContactVars& cVars,
        ibContactVars& tVars
    );

    vector getcLVec()
    {
        return cLVec_;
    }

    vector gettLVec()
    {
        return tLVec_;
    }

    vector getFNe();

    vector getFA();

    vector getFNd();

    vector getFt(scalar deltaT);

    Tuple2<forces,forces>& get_out_force()
    {
        return outForce_;
    }

    void setVMInfo(boundBox& bBox, scalar sub_volumeV);

    void setVMInfo(const virtual_meshInfo& vmInfo);

    std::shared_ptr<virtual_meshInfo>& getVMInfo();

    void clearOutForces()
    {
        outForce_.first().F = vector::zero;
        outForce_.first().T = vector::zero;
        outForce_.second().F = vector::zero;
        outForce_.second().T = vector::zero;
    }

    bool contact_resolved()
    {
        bool reduceBool = contact_resolved_;
        reduce(reduceBool,orOp<bool>());

        return reduceBool;
    }

    void setResolvedContact(bool resolved)
    {
        contact_resolved_ = resolved;
    }

    scalar& getLc()
    {
        return Lc_;
    }

    void syncData();
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
