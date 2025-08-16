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
    class for Verlet lists to optimize contact detection

SourceFiles
    verletList.C
\*---------------------------------------------------------------------------*/

#ifndef verletList_H
#define verletList_H

#include "dictionary.H"
#include "fvCFD.H"
#include "verletPoint.H"
#include "verletContact.H"
#include "immersed_body.H"

#include <list>
#include <map>
#include <unordered_set>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

struct hashFunction
{
  size_t operator()(const std::pair<label, label> &x) const
  {
    return x.first ^ x.second;
  }
};


class verletList
{
    private:

    typedef DynamicList<label> DynamicLabelList;
    typedef DynamicList<scalar> DynamicScalarList;
    typedef DynamicList<vector> DynamicVectorList;
    typedef DynamicList<point> DynamicPointList;

    typedef std::pair<label, label> cPair;
    typedef std::map<cPair, verletContact> cPairVcontactMap;

    typedef std::unordered_set<cPair, hashFunction> cPairHasSet;

    std::list<std::shared_ptr<verletBox>> verletBoxes_;

    List<std::list<std::shared_ptr<verletPoint>>> verletLists_;
  
    List<cPairVcontactMap> cntNeighList_;
 
    cPairHasSet posCntList_;

    void swapVerletPoints(
        std::shared_ptr<verletPoint> a,
        std::shared_ptr<verletPoint> b,
        label coord
    );

    void addCPairToCntNList(
        cPair cPair,
        label coord,
        std::shared_ptr<verletPoint> a,
        std::shared_ptr<verletPoint> b
    );

    public:

    // Constructors ------------------------------------------------
    verletList();

    // Destructors -------------------------------------------------
    ~verletList();

    void add_body_to_vList(immersed_body& ib);

    void removeBodyFromVList(immersed_body& ib);

    void initialSorting();

    void update(PtrList<immersed_body>& ibs);

    cPairHasSet::iterator begin()
    {
        return posCntList_.begin();
    }

    cPairHasSet::iterator end()
    {
        return posCntList_.end();
    }
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
