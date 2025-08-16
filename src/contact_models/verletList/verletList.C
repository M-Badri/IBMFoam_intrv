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
\*---------------------------------------------------------------------------*/
#include "verletList.H"

using namespace Foam;

//---------------------------------------------------------------------------//
verletList::verletList()
:
verletLists_(3),
cntNeighList_(3)
{}
verletList::~verletList()
{}
//---------------------------------------------------------------------------//
void verletList::swapVerletPoints(
    std::shared_ptr<verletPoint> a,
    std::shared_ptr<verletPoint> b,
    label coord
)
{
    if(a->getBodyId() != b->getBodyId())
    {
        if(!a->isMin() && b->isMin())
        {
            cPair newPair(
                min(a->getBodyId(), b->getBodyId()),
                max(a->getBodyId(), b->getBodyId())
            );
            addCPairToCntNList(newPair, coord, a, b);
        }
        else if(a->isMin() && !b->isMin())
        {
            cPair cPair(
                min(a->getBodyId(), b->getBodyId()),
                max(a->getBodyId(), b->getBodyId())
            );

            if(cntNeighList_[coord].find(cPair) != cntNeighList_[coord].end())
            {
                if (a->getBodyId() < b->getBodyId())
                {
                    cntNeighList_[coord][cPair].removeContactBox(a->getParentBox(), b->getParentBox());
                }
                else
                {
                    cntNeighList_[coord][cPair].removeContactBox(b->getParentBox(), a->getParentBox());
                }

                if (cntNeighList_[coord][cPair].getContactBoxes().size() == 0)
                {
                    cntNeighList_[coord].erase(cPair);
                }
            }

            if(posCntList_.find(cPair) != posCntList_.end() &&
                ((cntNeighList_[0].find(cPair) == cntNeighList_[0].end() && emptyDim != 0) ||
                (cntNeighList_[1].find(cPair) == cntNeighList_[1].end() && emptyDim != 1) ||
                (cntNeighList_[2].find(cPair) == cntNeighList_[2].end() && emptyDim != 2)))
            {
                posCntList_.erase(cPair);
            }
        }
    }
}
//---------------------------------------------------------------------------//
void verletList::addCPairToCntNList(
    cPair cPair,
    label coord,
    std::shared_ptr<verletPoint> a,
    std::shared_ptr<verletPoint> b
)
{
    if(cntNeighList_[coord].find(cPair) == cntNeighList_[coord].end())
    {
        verletContact newContact(cPair);
        if (a->getBodyId() < b->getBodyId())
        {
            newContact.addContactBox(a->getParentBox(), b->getParentBox());
        }
        else
        {
            newContact.addContactBox(b->getParentBox(), a->getParentBox());
        }

        cntNeighList_[coord].insert({cPair, newContact});
    }
    else
    {
        if (a->getBodyId() < b->getBodyId())
        {
            cntNeighList_[coord][cPair].addContactBox(a->getParentBox(), b->getParentBox());
        }
        else
        {
            cntNeighList_[coord][cPair].addContactBox(b->getParentBox(), a->getParentBox());
        }
    }

    if(a->getIsStatic() && b->getIsStatic())
    {
        return;
    }
    if(posCntList_.find(cPair) == posCntList_.end() &&
        (cntNeighList_[0].find(cPair) != cntNeighList_[0].end() || emptyDim == 0) &&
        (cntNeighList_[1].find(cPair) != cntNeighList_[1].end() || emptyDim == 1) &&
        (cntNeighList_[2].find(cPair) != cntNeighList_[2].end() || emptyDim == 2))
    {
        posCntList_.insert(cPair);
    }
}
//---------------------------------------------------------------------------//
void verletList::add_body_to_vList(immersed_body& ib)
{
    List<std::shared_ptr<boundBox>> bBoxes = ib.get_geom_model().getBBoxes();

    forAll(bBoxes, bBox)
    {
        verletBoxes_.push_back(verletBox::create(
            ib.getBodyId(),
            bBoxes[bBox],
            ib.getbodyOperation()==0
        ));

        verletBoxes_.back()->setVerletPoints();

        forAll(verletLists_, vListI)
        {
            if(vListI != emptyDim)
            {
                verletLists_[vListI].push_back(
                    verletBoxes_.back()->getMinPoint()
                );

                verletLists_[vListI].push_back(
                    verletBoxes_.back()->getMaxPoint()
                );
            }
        }
    }
}
//---------------------------------------------------------------------------//
void verletList::removeBodyFromVList(immersed_body& ib)
{
    forAll (verletLists_, coordI)
    {
        verletLists_[coordI].remove_if(
            [&ib](std::shared_ptr<verletPoint>& vPoint)
            {
                return vPoint->getBodyId() == ib.getBodyId();
            }
        );

        for (auto it = cntNeighList_[coordI].begin();
            it != cntNeighList_[coordI].end();)
        {
            if (it->first.first == ib.getBodyId()
                ||
                it->first.second == ib.getBodyId())
            {
                it = cntNeighList_[coordI].erase(it);
            }
            else
            {
                ++it;
            }
        }
    }

    verletBoxes_.remove_if(
        [&ib](std::shared_ptr<verletBox>& vBox)
        {
            return vBox->getBodyId() == ib.getBodyId();
        }
    );

    for (auto iter = posCntList_.begin(); iter != posCntList_.end();)
    {
        if (iter->first == ib.getBodyId()
            ||
            iter->second == ib.getBodyId())
        {
            iter = posCntList_.erase(iter);
        }
        else
        {
            ++iter;
        }
    }
}
//---------------------------------------------------------------------------//
void verletList::initialSorting()
{
    for (label coord = 0; coord < 3; ++coord)
    {
        verletLists_[coord].sort([&coord](
            std::shared_ptr<verletPoint>& vPoint1,
            std::shared_ptr<verletPoint>& vPoint2)
        {
            return vPoint1->getPoint()[coord] < vPoint2->getPoint()[coord];
        });

        std::list<std::shared_ptr<verletPoint>> openedIb;
        for (auto iter = verletLists_[coord].begin();
            iter != verletLists_[coord].end(); ++iter)
        {
            if((*iter)->isMin())
            {
                label curIb = (*iter)->getBodyId();
                for (auto oIbIter = openedIb.begin();
                    oIbIter != openedIb.end(); ++oIbIter)
                {
                    cPair newPair(
                        min(curIb, (*oIbIter)->getBodyId()), max(curIb, (*oIbIter)->getBodyId())
                    );
                    addCPairToCntNList(newPair, coord, *iter, *oIbIter);
                }

                openedIb.push_back(*iter);
            }
            else
            {
                label curIb = (*iter)->getBodyId();
                openedIb.remove_if(
                    [&curIb](std::shared_ptr<verletPoint>& vPoint)
                    {
                        return vPoint->getBodyId() == curIb;
                    }
                );
            }
        }
    }
}
//---------------------------------------------------------------------------//
void verletList::update(PtrList<immersed_body>& ibs)
{
    forAll(ibs, ibi)
    {
        ibs[ibi].get_geom_model().getBBoxes();
    }

    for (label coord = 0; coord < 3; ++coord)
    {
        if (!verletLists_[coord].empty())
        {
            auto it2 = verletLists_[coord].begin();
            auto it1 = it2++;

            while(true)
            {
                if((*it1)->getPoint()[coord]
                    > (*it2)->getPoint()[coord])
                {
                    swapVerletPoints(*it1, *it2, coord);
                    std::swap(*it1, *it2);

                    if (it1 != verletLists_[coord].begin())
                    {
                        it2 = it1--;
                        continue;
                    }
                }
                it1 = it2++;
                if (it2 == verletLists_[coord].end())
                {
                    break;
                }
            }
        }
    }
}
//---------------------------------------------------------------------------//
