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
#include "sub_contact.H"

using namespace Foam;

//---------------------------------------------------------------------------//
sub_contact::sub_contact()
{}

sub_contact::~sub_contact()
{}
//---------------------------------------------------------------------------//
void sub_contact::addSubVolume(std::shared_ptr<sub_volume> sV)
{
    if (sub_volumes_.size() == 0)
    {
        boundBox_ = *sV;
    }

    sub_volumes_.push_back(sV);
    tmp<pointField> points = boundBox_.points();
    points->append(sV->points());
    boundBox_ = boundBox(points,false);
    volume_ += sV->volume();
}
//---------------------------------------------------------------------------//
bool sub_contact::canCombine(sub_volume& sV)
{
    if (!boundBox_.overlaps(sV))
    {
        return false;
    }

    for (auto& sub_volume : sub_volumes_)
    {
        if (sub_volume->overlaps(sV))
        {
            return true;
        }
    }

    return false;
}
//---------------------------------------------------------------------------//
DynamicList<point> sub_contact::getEdgePoints() const
{
    DynamicList<point> edgePoints;
    for (auto& sub_volume : sub_volumes_)
    {
        if (sub_volume->isEdge())
        {
            edgePoints.append(sub_volume->midpoint());
        }
    }

    return edgePoints;
}
//---------------------------------------------------------------------------//
