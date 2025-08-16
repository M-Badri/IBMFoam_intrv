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
    prtContactInfo.C
\*---------------------------------------------------------------------------*/

#ifndef prtContactInfo_H
#define prtContactInfo_H

#include "dictionary.H"
#include "fvCFD.H"

#include "geom_model.H"
#include "ibContactVars.H"
#include "ibContactClass.H"
#include "prtSubContactInfo.H"
#include "physicalProperties.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{


class prtContactInfo
{
    private:

    typedef DynamicList<label> DynamicLabelList;
    typedef DynamicList<scalar> DynamicScalarList;
    typedef DynamicList<vector> DynamicVectorList;
    typedef DynamicList<point> DynamicPointList;

  
    protected:

    ibContactClass& cIbContactClass_;
    ibContactVars& cContactVars_;

    ibContactClass& tIbContactClass_;
    ibContactVars& tContactVars_;

    Tuple2<label,label> contactPair_;

    std::vector<std::shared_ptr<prtSubContactInfo>> contactList_;
    std::vector<std::shared_ptr<prtSubContactInfo>> newContactList_;

    physicalProperties physicalProperties_;

    std::shared_ptr<prtSubContactInfo> matchSubContact
    (
        boundBox& bbox,
        physicalProperties& physicalProperties,
        Tuple2<label,label>& contactPair
    );

    void limitBBox(boundBox& bbox);

    public:

    // Constructors ------------------------------------------------
    prtContactInfo
    (
        ibContactClass& cClass,
        ibContactVars& cVars,
        ibContactClass& tClass,
        ibContactVars& tVars
    );

    // Destructors -------------------------------------------------
    ~prtContactInfo();
    //-----------------------------Func-------------------------------//

    const Tuple2<label, label>& getCPair()
    {
        return contactPair_;
    }

    ibContactClass& getcClass()
    {
        return cIbContactClass_;
    }

    ibContactClass& gettClass()
    {
        return tIbContactClass_;
    }

    ibContactVars& getcVars()
    {
        return cContactVars_;
    }

    ibContactVars& gettVars()
    {
        return tContactVars_;
    }

    scalar getMu()
    {
        return physicalProperties_.aMu_;
    }

    void getContacts_Sphere();

    void getContacts_ArbShape
    (
        scalar cellV
    );

    void swapContactLists();

    bool contact_resolved();

    std::vector<std::shared_ptr<prtSubContactInfo>>& getPrtSCList()
    {
        return contactList_;
    }

    void syncContactList();

    void registerContactList(DynamicList<prtSubContactInfo*>& contactList);

    void clearData();

    void syncData();
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
