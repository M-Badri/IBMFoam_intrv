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

SourceFiles
    IBMFoam.C
\*---------------------------------------------------------------------------*/

#ifndef IBMFoam_H
#define IBMFoam_H

#include "dictionary.H"
#include "coordinateSystem.H"
#include "coordinateSystems.H"
#include "wordList.H"
#include "labelList.H"
#include "dimensionedScalar.H"
#include "dimensionedTensor.H"
#include "primitiveFieldsFwd.H"
#include "volFieldsFwd.H"
#include "fvMatricesFwd.H"
#include "triSurfaceMesh.H"
#include "immersed_body.H"

#include "add_model.H"
#include "add_modelOnce.H"
#include "add_modelOnceFromFile.H"
#include "add_modelOnceScatter.H"
#include "add_modelRepeatSamePosition.H"
#include "add_modelRepeatRandomPosition.H"
#include "add_modelDistribution.H"
#include "LS_interpol.H"
#include "line_interpol.H"
#include "verletList.H"

#include "clockTime.H"
namespace Foam
{

 class IBMFoam
 {

    private:
    typedef DynamicList<label>  DynamicLabelList;
    typedef DynamicList<scalar> DynamicScalarList;
    typedef DynamicList<vector> DynamicVectorList;
    typedef DynamicList<point>  DynamicPointList;

    const  fvMesh&          mesh_;
    IOdictionary            IBMFoamDict_;
    dictionary              IBMFoamInterpDict_;
    IOdictionary            trans_properties_;

    const wordList          body_names_;
    PtrList<immersed_body> imm_bodies_;

    autoPtr<ib_interpol> ib_interp_;
    PtrList<add_model> add_models_;

    List<labelList> cellPoints_;

    HashTable<autoPtr<prtContactInfo>,Tuple2<label, label>,Hash<Tuple2<label, label>>>
        prtcInfoTable_;

    scalar stepDEM_;

    label recomputeM0_;

    vector geometricD_;

    verletList verletList_;

    bool save_simulation_;
    bool rec_first_timeStep_ = false;
    word record_outDir_;

    public:

    // Constructors ----------------------------------------------------
    IBMFoam(const Foam::fvMesh& mesh);

    // Destructors -----------------------------------------------------
    ~IBMFoam();

    void initialize(volScalarField& body,volVectorField& U,volScalarField& refine_F,label recomputeM0, word runTime);

    void crt_bodies
    (
        volScalarField& body,
        volScalarField& refine_F
    );

    void pre_update_bodies
    (
        volScalarField& body,
        volVectorField& f
    );
    void post_update_bodies(volScalarField& body,volVectorField& f);

    void updateFSCoupling
    (
        volScalarField& body,
        volVectorField& f
    );

    void recrt_bodies(volScalarField& body,volScalarField& refine_F);
   
    void update_dem_(volScalarField& body,volScalarField& refine_F);

    prtContactInfo& getPrtcInfo(Tuple2<label,label> cPair);

    void interpolateIB
    (
        volVectorField& V,
        volVectorField& Vs,
        volScalarField& body
    );

    void writeBodiesInfo();

    void add_remove_bodies(volScalarField& body,volVectorField& U,volScalarField& refine_F);

    void restartSimulation
    (
        volScalarField& body,
        volScalarField& refine_F,
        word runTime
    );

    void preCalculateCellPoints();

    bool getRecordFirstTime()
    {
        return rec_first_timeStep_;
    }

    void setRecordFirstTime( bool setUp)
    {
        rec_first_timeStep_ = setUp;
    }

    void writeFirtsTimeBodiesInfo();
    
    void setSolverInfo();

 };
}

#endif
