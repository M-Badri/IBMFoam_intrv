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

#include "IBMFoam.H"
#include "polyMesh.H"
#include "fvCFD.H"
#include "fvMatrices.H"
#include "geometricOneField.H"

#include "interpolationCellPoint.H"
#include "interpolationCell.H"

#include "scalarMatrices.H"
#include "OFstream.H"
#include <iostream>
#include "def_externVars.H"
#include "parameters.H"

#define ORDER 2

using namespace Foam;
using namespace contactModel;

//---------------------------------------------------------------------------//
IBMFoam::IBMFoam(const Foam::fvMesh& mesh)
:
mesh_(mesh),
IBMFoamDict_
(
    IOobject
    (
        "IBMFoamDict",
        "constant",
        mesh_,
        IOobject::MUST_READ,
        IOobject::NO_WRITE
    )
),
trans_properties_
(
    IOobject
    (
        "transportProperties",
        "constant",
        mesh_,
        IOobject::MUST_READ,
        IOobject::NO_WRITE
    )
),
body_names_(IBMFoamDict_.lookup("body_names")),
prtcInfoTable_(0),
stepDEM_(readScalar(IBMFoamDict_.lookup("stepDEM"))),
save_simulation_(readBool(IBMFoamDict_.lookup("save_simulation")))
{
    materialProperties::matProps_insert(
        "None",
        material_info("None", 1, 1, 1, 1, 1)
    );

    if(IBMFoamDict_.found("recordFirstTimeStep"))
    {
        rec_first_timeStep_ = readBool(IBMFoamDict_.lookup("recordFirstTimeStep"));
    }

    if(IBMFoamDict_.found("nSolidsInDomain"))
    {
        solver_info::setNSolidsTreshnold(readLabel(IBMFoamDict_.lookup("nSolidsInDomain")));
    }

    dictionary dem_dict = IBMFoamDict_.subDict("DEM");
    dictionary materialsDic = dem_dict.subDict("materials");
    List<word> materials_names = materialsDic.toc();
    if(dem_dict.found("increasedDamping"))
    {
        contact_model_info::set_increased_damping(readBool(dem_dict.lookup("increasedDamping")));
    }
    else
    {
        contact_model_info::set_increased_damping(false);
    }

    forAll(materials_names, matI)
    {
        dictionary matI_dic = materialsDic.subDict(materials_names[matI]);
        scalar eps = readScalar(matI_dic.lookup("eps"));
        if(!contact_model_info::getIncreasedDamping())
        {
            eps = 0.906463027*eps + 0.093538298;//LinearRegression based on LIGGGHTS data testing
            eps = min(eps, 1.0);
        }

        materialProperties::matProps_insert(
            materials_names[matI],
            material_info(
                materials_names[matI],
                readScalar(matI_dic.lookup("Y")),
                readScalar(matI_dic.lookup("nu")),
                readScalar(matI_dic.lookup("mu")),
                readScalar(matI_dic.lookup("adh_coef")),
                eps
            )
        );
    }

    if(dem_dict.found("face_adh"))
    {
        dictionary interf_adh_dic = dem_dict.subDict("face_adh");
        List<word> interNames = interf_adh_dic.toc();
        forAll(interNames, interI)
        {
            dictionary interDicI = interf_adh_dic.subDict(interNames[interI]);
            wordList interMat = interDicI.lookup("materials");
            string interKey;
            if(interMat[0] < interMat[1])
            {
                interKey += interMat[0];
                interKey += "-";
                interKey += interMat[1];
            }
            else
            {
                interKey += interMat[1];
                interKey += "-";
                interKey += interMat[0];
            }

            inter_adhesion::inter_adhesion_insert(
                interKey,
                readScalar(interDicI.lookup("value"))
            );
        }
    }

    if(dem_dict.found("LcCoeff"))
    {
        contact_model_info::setLcCoeff(readScalar(dem_dict.lookup("LcCoeff")));
    }
    else
    {
        contact_model_info::setLcCoeff(4.0);
    }

    if(dem_dict.found("rotationModel"))
    {
        word rotModel = dem_dict.lookup("rotationModel");
        if(rotModel == "chen2012")
        {
            contact_model_info::set_rotation_model(0);
        }
        else if(rotModel == "mindlin1953")
        {
            contact_model_info::set_rotation_model(1);
        }
        else
        {
            Info << "Rotation Model not recognized, setting to default mindlin1953" << endl;
            contact_model_info::set_rotation_model(1);
        }    
    }
    else
    {
        Info << "Rotation Model not recognized, setting to default mindlin1953" << endl;
        contact_model_info::set_rotation_model(1);
    }
    

    Info <<" -- Coefficient for characteristic Lenght Lc is set to : "<< contact_model_info::getLcCoeff() << endl;

    dictionary patch_dict = dem_dict.subDict("coll_wall");
    List<word> patch_names = patch_dict.toc();
    forAll(patch_names, patchI)
    {
        word patchMaterial = patch_dict.subDict(patch_names[patchI]).lookup("material");
        vector patchNVec = patch_dict.subDict(patch_names[patchI]).lookup("nVec");
        vector pln_point = patch_dict.subDict(patch_names[patchI]).lookup("pln_point");

        wall_plane_info::wall_plane_info_insert(
            patch_names[patchI],
            patchNVec,
            pln_point
        );

        wall_mat_info::wall_mat_info_insert(
            patch_names[patchI],
            materialProperties::getMatProps()[patchMaterial]
        );
    }

    if(dem_dict.found("cyclicPatches"))
    {
        Info << "CyclicPatches Found " << endl;
        dictionary cyclicPatchDic = dem_dict.subDict("cyclicPatches");
        List<word> cyclic_patch_names = cyclicPatchDic.toc();
        forAll(cyclic_patch_names, patchI)
        {
            vector patchNVec = cyclicPatchDic.subDict(cyclic_patch_names[patchI]).lookup("nVec");
            vector pln_point = cyclicPatchDic.subDict(cyclic_patch_names[patchI]).lookup("pln_point");
            word neighbourPatch = cyclicPatchDic.subDict(cyclic_patch_names[patchI]).lookup("neighbourPatch");

            cyclic_plane_info::insert(
                cyclic_patch_names[patchI],
                patchNVec,
                pln_point,
                neighbourPatch
            );
        }
        Info << "Cyclic Patches  " <<  cyclic_patch_names <<endl;
    }

    if (IBMFoamDict_.found("geometricD"))
    {
        geometricD = IBMFoamDict_.lookup("geometricD");
    }
    else
    {
        geometricD = mesh_.geometricD();
    }

    forAll (geometricD, direction)
    {
        if (geometricD[direction] == -1)
        {
            case3D = false;
            emptyDir[direction] = 1;
            emptyDim = direction;
            break;
        }
    }

    if (IBMFoamDict_.isDict("virtual_mesh"))
    {
        dictionary vMDic = IBMFoamDict_.subDict("virtual_mesh");
        virtual_meshLevel::setVirtualMeshLevel(readScalar(vMDic.lookup("level")),readScalar(vMDic.lookup("charCellSize")));
        Info <<" -- Virt Mesh Decomposition Level is set to        : "<< virtual_meshLevel::getVirtualMeshLevel() << endl;
        Info <<" -- Virt Mesh char CellSize for boundary is set to  : "<< virtual_meshLevel::getCharCellSize() << endl;

    }
    else
    {
        virtual_meshLevel::setVirtualMeshLevel(1,1);
        Info <<" -- VirtMesh Decomposition Level is set to        : "<< virtual_meshLevel::getVirtualMeshLevel() << endl;
        Info <<" -- VirtMesh char CellSize for boundary is set to  : "<< virtual_meshLevel::getCharCellSize() << endl;

    }

    record_outDir_ = mesh_.time().rootPath() + "/" + mesh_.time().globalCaseName() + "/bodiesInfo";
}
//---------------------------------------------------------------------------//
IBMFoam::~IBMFoam()
{}
//---------------------------------------------------------------------------//
void IBMFoam::initialize
(
    volScalarField& body,
    volVectorField& U,
    volScalarField& refine_F,
    label recomputeM0,
    word runTime
)
{
    if(IBMFoamDict_.found("out_setting"))
    {
        dictionary outputDic = IBMFoamDict_.subDict("out_setting");
        bool basicOutput = readBool(outputDic.lookup("basic"));
        bool iBoutput = readBool(outputDic.lookup("iB"));
        bool DEMoutput = readBool(outputDic.lookup("DEM"));
        bool add_modelOutput = readBool(outputDic.lookup("add_model"));
        bool parallelDEMOutput = readBool(outputDic.lookup("parallelDEM"));
        InfoH.setOutput(
            basicOutput,
            iBoutput,
            DEMoutput,
            add_modelOutput,
            parallelDEMOutput
        );
    }

    preCalculateCellPoints();

    if(IBMFoamDict_.found("interpolationSchemes"))
    {
        IBMFoamInterpDict_ = IBMFoamDict_.subDict("interpolationSchemes");

        if(IBMFoamInterpDict_.found("method"))
        {
            word intMethod = IBMFoamInterpDict_.lookup("method");

            if(intMethod == "leastSquares")
            {
                dictionary lsCoeffsDict
                    = IBMFoamInterpDict_.subDict("leastSquaresCoeffs");
                ib_interp_.set(new LS_interpol(
                    mesh_,
                    readScalar(lsCoeffsDict.lookup("distFactor")),
                    readScalar(lsCoeffsDict.lookup("radiusFactor")),
                    readScalar(lsCoeffsDict.lookup("angleFactor")),
                    readScalar(lsCoeffsDict.lookup("maxCCRows"))
                ));
            }
            else if(intMethod == "line")
            {
                ib_interp_.set(new line_interpol(IBMFoamInterpDict_));
            }
        }
    }

    bool startTime0(runTime == "0");

    // initialize add_models
    add_models_.setSize(body_names_.size());
    imm_bodies_.setSize(0);                                         
    refine_F *= 0;
    recomputeM0_ = recomputeM0;

    if(!startTime0)
    {
        if(!isDir(record_outDir_))
            mkDir(record_outDir_);
        else
        {
            fileNameList entries(readDir(record_outDir_,fileType::directory)); 
            scalar runTimeS(stod(runTime));
            forAll(entries,entry)
            {
                scalar dirTime(stod(entries[entry].name()));
                if(dirTime > runTimeS)
                {
                    word pathI(record_outDir_ + "/" + entries[entry]);
                    rmDir(pathI);
                }
            }
        }

        restartSimulation(body, refine_F, runTime);
    }
    else
    {
        if(!isDir(record_outDir_))
            mkDir(record_outDir_);
        else
        {
            rmDir(record_outDir_);
            mkDir(record_outDir_);
        }
    }

    #include "initialize_add_models.H"

    forAll (add_models_,modelI)
    {
        word body_name(body_names_[modelI]);
        InfoH << basic_Info << "Creating immersed body based on: " << body_name << endl;

        label max_additions(1000);
        label c_Addition(0);

        while (add_models_[modelI].shouldAddBody(body) and c_Addition < max_additions and imm_bodies_.size() < solver_info::getNSolidsTreshnold())
        {
            InfoH << add_model_Info << "add_model invoked action, trying to add new body" << endl;
            std::shared_ptr<geom_model> b_geomModel(add_models_[modelI].addBody(body, imm_bodies_));
            c_Addition++;
            if (add_models_[modelI].getBodyAdded())
            {
                label newIBSize(imm_bodies_.size()+1);
                label add_IB_pos(newIBSize - 1);
                imm_bodies_.setSize(newIBSize);

                InfoH << add_model_Info << "Trying to set imm_bodies" << endl;
                imm_bodies_.set
                (
                    add_IB_pos,
                    new immersed_body
                    (
                        body_name,
                        mesh_,
                        IBMFoamDict_,
                        trans_properties_,
                        add_IB_pos,
                        recomputeM0_,
                        b_geomModel,
                        ib_interp_,
                        cellPoints_
                    )
                );
                imm_bodies_[add_IB_pos].create_immersed_body(body,refine_F);
                imm_bodies_[add_IB_pos].compute_body_charPars();
                if (imm_bodies_[add_IB_pos].get_start_synced())
                {
                    imm_bodies_[add_IB_pos].initSyncWithFlow(U);
                }
                verletList_.add_body_to_vList(imm_bodies_[add_IB_pos]);
                InfoH << add_model_Info << "Body based on: " << body_name << " successfully added" << endl;
                c_Addition = 0;
            }
            else
            {
                InfoH << add_model_Info << "Body based on: "
                    << body_name << " should have been added but was not "
                    << "(probably overlap with an already existing body)"
                    << endl;
            }
        }
    }

    verletList_.initialSorting();
}
//---------------------------------------------------------------------------//
void IBMFoam::crt_bodies(volScalarField& body,volScalarField& refine_F)
{
    forAll (imm_bodies_,body_id)
    {
        if (imm_bodies_[body_id].get_is_active())
        {
            imm_bodies_[body_id].postContactUpdateBodyField(body,refine_F);
        }
    }

    DynamicList<scalar> particle_masses;
    DynamicList<label> particle_cells;
    DynamicList<symmTensor> particle_inertia_tensors;

    forAll (imm_bodies_,body_id)
    {
        if (imm_bodies_[body_id].get_is_active())
        {
            imm_bodies_[body_id].syncImmersedBodyParralell1(body,refine_F);
            if (imm_bodies_[body_id].get_geom_model().isCluster())
            {
                cluster_body& cBody = dynamic_cast<cluster_body&>(imm_bodies_[body_id].get_geom_model());
                std::vector<std::shared_ptr<geom_model>>& cBodies = cBody.get_cluster_bodies();
                for (auto& cB : cBodies)
                {
                    particle_masses.append(cB->getM());
                    particle_cells.append(cB->getNCells());
                    particle_inertia_tensors.append(cB->getI());
                }
            }
            else
            {
                particle_masses.append(imm_bodies_[body_id].get_geom_model().getM());
                particle_cells.append(imm_bodies_[body_id].get_geom_model().getNCells());
                particle_inertia_tensors.append(imm_bodies_[body_id].get_geom_model().getI());
            }
        }
    }
    reduce(particle_masses,sumOp<List<scalar>>());
    reduce(particle_cells,sumOp<List<label>>());
    reduce(particle_inertia_tensors,sumOp<List<symmTensor>>());

    label bodyIndex(0);
    forAll (imm_bodies_,body_id)
    {
        if (imm_bodies_[body_id].get_is_active())
        {
            if (imm_bodies_[body_id].get_geom_model().isCluster())
            {
                cluster_body& cBody = dynamic_cast<cluster_body&>(imm_bodies_[body_id].get_geom_model());
                std::vector<std::shared_ptr<geom_model>>& cBodies = cBody.get_cluster_bodies();
                for (auto& cB : cBodies)
                {
                    cB->setM(particle_masses[bodyIndex]);
                    cB->setNCells(particle_cells[bodyIndex]);
                    cB->setI(particle_inertia_tensors[bodyIndex]);
                    bodyIndex++;
                }
                cBody.setMassAndInertia();
            }
            else
            {
                imm_bodies_[body_id].get_geom_model().setM(particle_masses[bodyIndex]);
                imm_bodies_[body_id].get_geom_model().setNCells(particle_cells[bodyIndex]);
                imm_bodies_[body_id].get_geom_model().setI(particle_inertia_tensors[bodyIndex]);
                bodyIndex++;
            }

            imm_bodies_[body_id].syncImmersedBodyParralell2(body,refine_F);
            imm_bodies_[body_id].check_if_in_domain(body);
            imm_bodies_[body_id].update_old_movement_vars();
        }
    }

    forAll (imm_bodies_,body_id)
    {
        if (imm_bodies_[body_id].get_is_active())
        {
            imm_bodies_[body_id].chceckBodyOp();
        }
    }
}
//---------------------------------------------------------------------------//
void IBMFoam::pre_update_bodies
(
    volScalarField& body,
    volVectorField& f
)
{
    forAll (imm_bodies_,body_id)
    {
        if (imm_bodies_[body_id].get_is_active())
        {
            imm_bodies_[body_id].inContact_with_static(false);

            imm_bodies_[body_id].update_old_movement_vars();
            imm_bodies_[body_id].printStats();
        }
    }
}
//---------------------------------------------------------------------------//
void IBMFoam::post_update_bodies
(
    volScalarField& body,
    volVectorField& f
)
{
    forAll (imm_bodies_,body_id)
    {
        if (imm_bodies_[body_id].get_is_active())
        {
            imm_bodies_[body_id].clearIntpInfo();
            imm_bodies_[body_id].post_pimple_update_immersed_body(body,f);
        }
    }
}
//---------------------------------------------------------------------------//
void IBMFoam::recrt_bodies
(
    volScalarField& body,
    volScalarField& refine_F
)
{
    refine_F *= 0;
    preCalculateCellPoints();
    forAll (add_models_,modelI)
    {
        add_models_[modelI].recreateBoundBox();
    }
    forAll (imm_bodies_,body_id)
    {
        if (imm_bodies_[body_id].get_is_active())
        {
            imm_bodies_[body_id].recreateBodyField(body,refine_F);
        }
    }
    forAll (imm_bodies_,body_id)
    {
        if (imm_bodies_[body_id].get_is_active())
        {
            imm_bodies_[body_id].syncCreateImmersedBody(body,refine_F);
            imm_bodies_[body_id].check_if_in_domain(body);
            if(imm_bodies_[body_id].getrecomputeM0() > 0)
            {
                imm_bodies_[body_id].compute_body_charPars();
                imm_bodies_[body_id].recomputedM0();
            }
            InfoH << iB_Info << "-- body "
                << imm_bodies_[body_id].getBodyId() << " Re-created" << endl;
        }
    }
}
//---------------------------------------------------------------------------//
void IBMFoam::interpolateIB( volVectorField & V
                              ,volVectorField & Vs
                              ,volScalarField & body)
{
    if(ib_interp_.valid())
    {
        ib_interp_->resetInterpolator(V);
    }
    Vs = V;

    forAll (imm_bodies_,body_id)
    {
        if (imm_bodies_[body_id].get_is_active())
        {
            imm_bodies_[body_id].update_vector_field(Vs, V.name(),body);

            if(ib_interp_.valid())
            {
                ib_interp_->ib_interpolate
                (
                    imm_bodies_[body_id].getIntpInfo(),
                    Vs,
                    imm_bodies_[body_id].getUatIbPoints(),
                    mesh_
                );
            }
        }
    }
}
//---------------------------------------------------------------------------//
void IBMFoam::writeBodiesInfo()
{
    if(!save_simulation_)
        return;

    word curOutDir(record_outDir_ + "/" + mesh_.time().timeName());


    mkDir(curOutDir);
    mkDir(curOutDir +"/stlFiles");
    DynamicLabelList active_ib;
    forAll (imm_bodies_,body_id)
    {
        if (imm_bodies_[body_id].get_is_active())
        {
            active_ib.append(body_id);
        }
    }
    wordList body_names;
    scalar listZize(active_ib.size());
    label bodies_per_proc = ceil(listZize/Pstream::nProcs());
    InfoH << basic_Info << "Active IB listZize      : " << listZize<< endl;
    InfoH << basic_Info << "bodies_per_proc : " << bodies_per_proc<< endl;

    for(int assignProc = Pstream::myProcNo()*bodies_per_proc; assignProc < min((Pstream::myProcNo()+1)*bodies_per_proc,active_ib.size()); assignProc++)
    {
        const label body_id(active_ib[assignProc]);
        word path(curOutDir + "/body" + std::to_string(imm_bodies_[body_id].getBodyId()) +".info");
        OFstream ofStream(path);
        IOobject outClass
            (
                path,
                mesh_,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            );
        IOdictionary outDict(outClass);

        outDict.writeHeader(ofStream);
        imm_bodies_[body_id].record_body_info(outDict,curOutDir);
        outDict.writeData(ofStream);
    }

}
//---------------------------------------------------------------------------//
void IBMFoam::update_dem_(volScalarField& body,volScalarField& refine_F)
{
    if (cyclic_plane_info::getCyclicPlaneInfo().size() > 0)
    {
        forAll (imm_bodies_,body_id)
        {
            if (!imm_bodies_[body_id].get_geom_model().isCluster())
            {
                vector transVec = vector::zero;

                if (detectCyclicContact(
                    imm_bodies_[body_id].get_wall_cnt_info(),
                    transVec
                ))
                {
                    verletList_.removeBodyFromVList(imm_bodies_[body_id]);

                    scalar thrSurf(readScalar(IBMFoamDict_.lookup("surface_threshold")));
                    std::shared_ptr<periodicBody> newPeriodicBody
                        = std::make_shared<periodicBody>(mesh_, thrSurf);

                    newPeriodicBody->setRhoS(imm_bodies_[body_id].get_geom_model().getRhoS());
                    std::shared_ptr<geom_model> iBcopy(imm_bodies_[body_id].get_geom_model().getCopy());
                    iBcopy->bodyMovePoints(transVec);
                    newPeriodicBody->add_body_to_cluster(imm_bodies_[body_id].get_geom_modelPtr());
                    newPeriodicBody->add_body_to_cluster(iBcopy);
                    imm_bodies_[body_id].get_geom_modelPtr() = newPeriodicBody;

                    verletList_.add_body_to_vList(imm_bodies_[body_id]);
                    Info << "Periodic body created for body " << body_id << endl;
                }
            }
            else
            {
                periodicBody& cBody = dynamic_cast<periodicBody&>(imm_bodies_[body_id].get_geom_model());

                if(cBody.shouldBeUnclustered())
                {
                    verletList_.removeBodyFromVList(imm_bodies_[body_id]);

                    imm_bodies_[body_id].get_geom_modelPtr() = cBody.getRemGeomModel();

                    verletList_.add_body_to_vList(imm_bodies_[body_id]);
                    Info << "Periodic body unclustered for body " << body_id << endl;
                }
            }
        }
    }

    scalar deltaTime(mesh_.time().deltaT().value());
    scalar pos(0.0);
    scalar step(stepDEM_);
    List<DynamicList<pointField>> bodies_position_list(Pstream::nProcs());
    HashTable <label,Tuple2<label, label>,Hash<Tuple2<label, label>>> sync_out_force_key_table;
    HashTable <label,Tuple2<label, label>,Hash<Tuple2<label, label>>> contact_resolvedKeyTable;
    HashTable <label,label,Hash<label>> wall_contact_ibTable;
    while( pos < 1)
    {
        bodies_position_list[Pstream::myProcNo()].clear();

        InfoH << DEM_Info << " Start DEM pos: " << pos
            << " DEM step: " << step << endl;

        InfoH << basic_Info << " DEM - CFD Time: "
            << mesh_.time().value() + deltaTime*pos << endl;

        forAll (imm_bodies_,ib)
        {
            imm_bodies_[ib].updateMovement(deltaTime*step*0.5);

            if(Pstream::myProcNo() == 0 )
            {
                imm_bodies_[ib].moveImmersedBody(deltaTime*step);
                if(imm_bodies_[ib].get_geom_model().getcType() != cluster)
                {
                    bodies_position_list[Pstream::myProcNo()].append(imm_bodies_[ib].get_geom_model().getBodyPoints());
                }
                else
                {
                    cluster_body& cBody = dynamic_cast<cluster_body&>(imm_bodies_[ib].get_geom_model());
                    std::vector<std::shared_ptr<geom_model>>& cBodies = cBody.get_cluster_bodies();
                    for (auto& cB : cBodies)
                    {
                        bodies_position_list[Pstream::myProcNo()].append(cB->getBodyPoints());
                    }
                }
            }
        }

        Pstream::gatherList(bodies_position_list,0);
        Pstream::scatterList(bodies_position_list,0);

        label bodyIndex(0);
        forAll (imm_bodies_,ib)
        {
            if(imm_bodies_[ib].get_geom_model().getcType() != cluster)
            {
                imm_bodies_[ib].get_geom_model().setBodyPosition(bodies_position_list[0][bodyIndex++]);
            }
            else
            {
                cluster_body& cBody = dynamic_cast<cluster_body&>(imm_bodies_[ib].get_geom_model());
                std::vector<std::shared_ptr<geom_model>>& cBodies = cBody.get_cluster_bodies();
                for (auto& cB : cBodies)
                {
                    cB->setBodyPosition(bodies_position_list[0][bodyIndex++]);
                }
            }
        }

        bodies_position_list[Pstream::myProcNo()].clear();

        verletList_.update(imm_bodies_);

        DynamicLabelList wall_contact_ib;
        wall_contact_ibTable.clear();
        forAll (imm_bodies_,body_id)
        {
            immersed_body& cIb(imm_bodies_[body_id]);
            if (cIb.get_is_active())
            {
                cIb.reset_contact_forces();

                if(cIb.getbodyOperation() != 0)
                {
                    if(detectWallContact
                    (
                        mesh_,
                        cIb.get_ib_contact_class(),
                        cIb.get_wall_cnt_info()
                    ))
                    {
                        cIb.get_ib_contact_class().setWallContact(true);
                        cIb.get_ib_contact_class().inContact_with_static(true);
                        wall_contact_ib.append(body_id);
                        wall_contact_ibTable.insert(body_id,wall_contact_ib.size()-1);
                    }
                }
            }
        }
        List<bool> wall_contact_resolved_list(wall_contact_ib.size(),false);

        if(wall_contact_ib.size() > 0)
        {
            label wallContactPerProc(ceil(double(wall_contact_ib.size())/Pstream::nProcs()));
            if( wall_contact_ib.size() <= Pstream::nProcs())
            {
                wallContactPerProc = 1;
            }
            for(int assignProc = Pstream::myProcNo()*wallContactPerProc; assignProc < min((Pstream::myProcNo()+1)*wallContactPerProc,wall_contact_ib.size()); assignProc++)
            {
                immersed_body& cIb(imm_bodies_[wall_contact_ib[assignProc]]);
                if(cIb.get_geom_model().getcType() != sphere && cIb.get_geom_model().getcType() != cluster)
                {
                    cIb.get_wall_cnt_info().findContactAreas();
                }

                DynamicList<wallSubContactInfo*> wall_contact_list;
                cIb.get_wall_cnt_info().registerSubContactList(wall_contact_list);
                List<bool> wallcRList(wall_contact_list.size(),false);

                forAll(wall_contact_list,sC)
                {
                    wallSubContactInfo* sCW = wall_contact_list[sC];
                    bool resolved(solveWallContact(
                        mesh_,
                        cIb.get_wall_cnt_info(),
                        deltaTime*step,
                        *sCW
                        ));
                    sCW->setResolvedContact(resolved);
                    wall_contact_resolved_list[assignProc] += resolved;
                    wallcRList[sC] = resolved;
                }
            }

            reduce(wall_contact_resolved_list,sumOp<List<bool>>());

            List<vector> iBody_OutForce_list(wall_contact_ib.size(),vector::zero);
            List<vector> iBody_OutTorque_list(wall_contact_ib.size(),vector::zero);

            forAll (wall_contact_ib,iB)
            {
                immersed_body& cIb(imm_bodies_[wall_contact_ib[iB]]);
                if(wall_contact_ibTable.found(cIb.getBodyId()))
                {
                    label cKey(wall_contact_ibTable[cIb.getBodyId()]);
                    if(wall_contact_resolved_list[cKey])
                    {
                        std::vector<std::shared_ptr<wallSubContactInfo>>& subCList
                            = cIb.get_wall_cnt_info().getWallSCList();

                        for(auto sCW : subCList)
                        {
                            iBody_OutForce_list[cKey] += sCW->get_out_force().F;
                            iBody_OutTorque_list[cKey] += sCW->get_out_force().T;
                        }
                    }
                }
            }
            reduce(iBody_OutForce_list,sumOp<List<vector>>());
            reduce(iBody_OutTorque_list,sumOp<List<vector>>());

            forAll (wall_contact_ib,iB)
            {
                immersed_body& cIb(imm_bodies_[wall_contact_ib[iB]]);
                forces cF;
                cF.F = iBody_OutForce_list[iB];
                cF.T = iBody_OutTorque_list[iB];

                cIb.updateContactForces(cF);
                cIb.get_wall_cnt_info().clearOldContact();
            }
        }

        wall_contact_ib.clear();

        DynamicList<prtSubContactInfo*> contactList;
        label vListSize(0);
        for (auto it = verletList_.begin(); it != verletList_.end(); ++it)
        {
            const Tuple2<label, label> cPair = Tuple2<label, label>(it->first, it->second);

            label cInd(cPair.first());
            bool cStatic(imm_bodies_[cInd].getbodyOperation() == 0);

            label tInd(cPair.second());
            bool tStatic(imm_bodies_[tInd].getbodyOperation() == 0);

            if((imm_bodies_[cInd].get_is_active() && imm_bodies_[tInd].get_is_active())
                &&
                !(cStatic && tStatic)
            )
            {
                if(cStatic)
                    imm_bodies_[tInd].inContact_with_static(true);

                if(tStatic)
                    imm_bodies_[cInd].inContact_with_static(true);

                prtContactInfo& prtcInfo(getPrtcInfo(
                    cPair)
                );

                prtcInfo.clearData();
                getContacts(
                    mesh_,
                    prtcInfo
                );

                prtcInfo.registerContactList(contactList);
            }
            vListSize++;
        }

        List<bool> contact_resolved(contactList.size(),false);
        List<label> contact_resolvedcKey(contactList.size(),0);
        List<label> contact_resolvedtKey(contactList.size(),0);
        bool syncedData(true);
        reduce(syncedData, orOp<bool>());

        if(contactList.size() > 0 )
        {
            label contactPerProc(ceil(double(contactList.size())/Pstream::nProcs()));
            if( contactList.size() <= Pstream::nProcs())
            {
                contactPerProc = 1;
            }

            for(int assignProc = Pstream::myProcNo()*contactPerProc; assignProc < min((Pstream::myProcNo()+1)*contactPerProc,contactList.size()); assignProc++)
            {
                prtSubContactInfo* sCI = contactList[assignProc];
                const Tuple2<label, label>& cPair = sCI->getCPair();

                contact_resolvedcKey[assignProc] = cPair.first();
                contact_resolvedtKey[assignProc] = cPair.second();

                ibContactClass& cClass(imm_bodies_[cPair.first()].get_ib_contact_class());
                ibContactClass& tClass(imm_bodies_[cPair.second()].get_ib_contact_class());

                if(detectPrtPrtContact(mesh_,cClass,tClass,*sCI))
                {
                    prtContactInfo& prtcInfo(getPrtcInfo(cPair));

                    bool resolved(solvePrtContact(mesh_, prtcInfo, *sCI, deltaTime*step));
                    sCI->setResolvedContact(resolved);

                    contact_resolved[assignProc] += resolved;
                }
            }
        }

        reduce(contact_resolved,sumOp<List<bool>>());
        reduce(contact_resolvedcKey,sumOp<List<label>>());
        reduce(contact_resolvedtKey,sumOp<List<label>>());
        contact_resolvedKeyTable.clear();
        forAll(contact_resolvedcKey,cKey)
        {
            contact_resolvedKeyTable.insert(Tuple2<label, label>(contact_resolvedcKey[cKey],contact_resolvedtKey[cKey]),cKey);
        }
        List<vector> cBodyOutForceList(vListSize,vector::zero);
        List<vector> cBodyOutTorqueList(vListSize,vector::zero);
        List<vector> tBodyOutForceList(vListSize,vector::zero);
        List<vector> tBodyOutTorqueList(vListSize,vector::zero);
        sync_out_force_key_table.clear();

        label nIter(0);
        for (auto it = verletList_.begin(); it != verletList_.end(); ++it)
        {
            const Tuple2<label, label> cPair = Tuple2<label, label>(it->first, it->second);

            prtContactInfo& prtcInfo(getPrtcInfo(cPair));

            if(contact_resolvedKeyTable.found(cPair))
            {
                label nSubContact(0);
                std::vector<std::shared_ptr<prtSubContactInfo>>& subCList
                    = prtcInfo.getPrtSCList();
                for(auto sC : subCList)
                {
                    nSubContact++;
                    cBodyOutForceList[nIter] += sC->get_out_force().first().F;
                    cBodyOutTorqueList[nIter] += sC->get_out_force().first().T;
                    tBodyOutForceList[nIter] += sC->get_out_force().second().F;
                    tBodyOutTorqueList[nIter] += sC->get_out_force().second().T;
                }
            }
            sync_out_force_key_table.insert(cPair,nIter);
            nIter++;
        }

        reduce(cBodyOutForceList,sumOp<List<vector>>());
        reduce(cBodyOutTorqueList,sumOp<List<vector>>());
        reduce(tBodyOutForceList,sumOp<List<vector>>());
        reduce(tBodyOutTorqueList,sumOp<List<vector>>());

        label nvListIter(0);

        for (auto it = verletList_.begin(); it != verletList_.end(); ++it)
        {
            const Tuple2<label, label> cPair = Tuple2<label, label>(it->first, it->second);
            label cInd(cPair.first());
            label tInd(cPair.second());

            if(!contact_resolvedKeyTable.found(cPair))
            {
                if(prtcInfoTable_.found(cPair))
                {
                    prtcInfoTable_.erase(cPair);
                    continue;
                }
            }
            else if(contact_resolved[contact_resolvedKeyTable[cPair]])
            {
                if(!sync_out_force_key_table.found(cPair))
                {
                    Pout <<" -- cPair  "<<cInd << " - "<<tInd << " not found in sync_out_force_key_table" << endl;
                    continue;
                }

                nvListIter = sync_out_force_key_table[cPair];
                if(nvListIter > cBodyOutForceList.size())
                {
                    Pout <<" -- cPair  "<<cInd << " - "<<tInd << " nvListIter > bodiesOutForceList[Pstream::myProcNo()].size()" << endl;
                    continue;
                }

                vector F1 = vector::zero;
                vector T1 = vector::zero;
                vector F2 = vector::zero;
                vector T2 = vector::zero;

                F1 += cBodyOutForceList[nvListIter];
                T1 += cBodyOutTorqueList[nvListIter];
                F2 += tBodyOutForceList[nvListIter];
                T2 += tBodyOutTorqueList[nvListIter];

                forces cF;
                cF.F = F1;
                cF.T = T1;
                forces tF;
                tF.F = F2;
                tF.T = T2;

                imm_bodies_[cInd].updateContactForces(cF);
                imm_bodies_[tInd].updateContactForces(tF);
            }
            else
            {
                if(prtcInfoTable_.found(cPair))
                {
                    prtcInfoTable_.erase(cPair);
                    continue;
                }
            }
        }

        forAll (imm_bodies_,ib)
        {
            imm_bodies_[ib].updateMovement(deltaTime*step*0.5);
            imm_bodies_[ib].printBodyInfo();

        }

        pos += step;

        if (pos + step + SMALL >= 1)
            step = 1 - pos;
    }
}
//---------------------------------------------------------------------------//
prtContactInfo& IBMFoam::getPrtcInfo(Tuple2<label,label> cPair)
{
    if(!prtcInfoTable_.found(cPair))
    {
        prtcInfoTable_.insert(cPair, autoPtr<prtContactInfo>( new prtContactInfo(
            imm_bodies_[cPair.first()].get_ib_contact_class(),
            imm_bodies_[cPair.first()].getContactVars(),
            imm_bodies_[cPair.second()].get_ib_contact_class(),
            imm_bodies_[cPair.second()].getContactVars()
        )));
    }

    return prtcInfoTable_[cPair]();
}
//---------------------------------------------------------------------------//

void IBMFoam::add_remove_bodies
(
    volScalarField& body,
    volVectorField& U,
    volScalarField& refine_F
)
{
    forAll (add_models_,modelI)
    {
        word body_name(body_names_[modelI]);

        label max_additions(50);
        label c_Addition(0);

        while (add_models_[modelI].shouldAddBody(body) and c_Addition < max_additions)
        {
            InfoH << add_model_Info << "add_model invoked action, trying to add new body" << endl;
            std::shared_ptr<geom_model> b_geomModel(add_models_[modelI].addBody(body, imm_bodies_));

            c_Addition++;

            if (add_models_[modelI].getBodyAdded())
            {
                InfoH << add_model_Info << "STL file correctly generated, registering the new body" << endl;
                label newIBSize(imm_bodies_.size()+1);
                label add_IB_pos(newIBSize - 1);
                imm_bodies_.setSize(newIBSize);
                imm_bodies_.set
                (
                    add_IB_pos,
                    new immersed_body
                    (
                        body_name,
                        mesh_,
                        IBMFoamDict_,
                        trans_properties_,
                        add_IB_pos,
                        recomputeM0_,
                        b_geomModel,
                        ib_interp_,
                        cellPoints_
                    )
                );
                immersed_body& nBody(imm_bodies_[add_IB_pos]);
                nBody.create_immersed_body(body,refine_F);
                nBody.compute_body_charPars();
                if (nBody.get_start_synced())
                {
                    nBody.initSyncWithFlow(U);
                }
                verletList_.add_body_to_vList(nBody);

                InfoH << add_model_Info
                    << "new body included into the simulation" << endl;
                c_Addition = 0;
            }
            else
            {
                InfoH << add_model_Info
                    << "new body should have been added but was not "
                    << "(probably overlap with an existing body)"
                    << endl;
            }
        }
    }
}
//---------------------------------------------------------------------------//
void IBMFoam::updateFSCoupling
(
    volScalarField& body,
    volVectorField& f
)
{
    forAll (imm_bodies_,body_id)
    {
        if (imm_bodies_[body_id].get_is_active())
        {
            imm_bodies_[body_id].pimpleUpdate(body,f);
        }
    }
}
//---------------------------------------------------------------------------//
void IBMFoam::restartSimulation
(
    volScalarField& body,
    volScalarField& refine_F,
    word runTime
)
{
    word timePath(record_outDir_+"/"+runTime);
    fileNameList files(readDir(timePath));
    scalar thrSurf(readScalar(IBMFoamDict_.lookup("surface_threshold")));

    forAll(files,f)
    {
        IOdictionary body_dict
        (
            IOobject
            (
                timePath + "/" + files[f],
                mesh_,
                IOobject::MUST_READ,
                IOobject::NO_WRITE
            )
        );

        word body_id(std::to_string(readLabel(body_dict.lookup("body_id"))));
        word body_name(body_dict.lookup("body_name"));
        vector Vel(body_dict.lookup("Vel"));
        scalar omega(readScalar(body_dict.lookup("omega")));
        vector Axis(body_dict.lookup("Axis"));
        bool isStatic(readBool(body_dict.lookup("static")));
        label timeStepsInContWStatic(readLabel(body_dict.lookup("timeStepsInContWStatic")));

        std::shared_ptr<geom_model> b_geomModel;
        word b_geom;
        if (IBMFoamDict_.subDict(body_name).found("b_geom"))
        {
            word input = IBMFoamDict_.subDict(body_name).lookup("b_geom");
            b_geom = input;
            InfoH << iB_Info << "Found b_geom for "
                << body_name << ", the body is: " << b_geom << endl;
        }
        else
        {
            b_geom = "convex";
            InfoH << iB_Info << "Did not find b_geom for "
                << body_name << ", using b_geom: " << b_geom << endl;
        }

        if(b_geom == "convex")
        {
            word stlPath(timePath + "/stlFiles/"+body_id+".stl");
            b_geomModel = std::make_shared<convex_body>(mesh_,stlPath,thrSurf);
        }
        else if(b_geom == "nonConvex")
        {
            word stlPath(timePath + "/stlFiles/"+body_id+".stl");
            b_geomModel = std::make_shared<nonConvex_body>(mesh_,stlPath,thrSurf);
        }
        else if(b_geom == "sphere")
        {
            vector startPosition = body_dict.subDict("sphere").lookup("position");
            scalar radius = readScalar(body_dict.subDict("sphere").lookup("radius"));

            b_geomModel = std::make_shared<sphere_body>(mesh_,startPosition,radius,thrSurf);
        }
        else
        {
            word stlPath(timePath + "/stlFiles/"+body_id+".stl");
            InfoH << iB_Info << "b_geom: " << b_geom
                << " not supported, using b_geom nonConvex" << endl;
            b_geom = "nonConvex";
            b_geomModel = std::make_shared<nonConvex_body>(mesh_,stlPath,thrSurf);
        }

        label newIBSize(imm_bodies_.size()+1);
        label add_IB_pos(newIBSize - 1);
        imm_bodies_.setSize(newIBSize);

        InfoH << iB_Info << "Restarting body: " << body_id << " as "
            << add_IB_pos << " body_name: " << body_name << endl;
        imm_bodies_.set
        (
            add_IB_pos,
            new immersed_body
            (
                body_name,
                mesh_,
                IBMFoamDict_,
                trans_properties_,
                add_IB_pos,
                recomputeM0_,
                b_geomModel,
                ib_interp_,
                cellPoints_
            )
        );

        imm_bodies_[add_IB_pos].create_immersed_body(body,refine_F);
        imm_bodies_[add_IB_pos].compute_body_charPars();
        imm_bodies_[add_IB_pos].setRestartSim(Vel,omega,Axis,isStatic,timeStepsInContWStatic);
        verletList_.add_body_to_vList(imm_bodies_[add_IB_pos]);
    }
}
//---------------------------------------------------------------------------//
void IBMFoam::preCalculateCellPoints()
{
    cellPoints_.clear();
    cellPoints_.setSize(mesh_.nCells());
    forAll(mesh_.C(), cellI)
    {
        cellPoints_[cellI] = mesh_.cellPoints()[cellI];
    }

    forAll (imm_bodies_,body_id)
    {
        imm_bodies_[body_id].get_geom_model().resetHashTable();
    }
}
//---------------------------------------------------------------------------//
void IBMFoam::writeFirtsTimeBodiesInfo()
{
    word curOutDir(record_outDir_ + "/" + mesh_.time().timeName());
    bool checkExistance(false);
    if(!save_simulation_ || isDir(curOutDir))
        return;
    reduce(checkExistance,orOp<bool>());
    if(Pstream::myProcNo() == 0)
    {
        mkDir(curOutDir);
        mkDir(curOutDir +"/stlFiles");
    }
    reduce(checkExistance,orOp<bool>());

    DynamicLabelList active_ib;
    forAll (imm_bodies_,body_id)
    {
        if (imm_bodies_[body_id].get_is_active())
        {
            active_ib.append(body_id);
        }
    }

    wordList body_names;
    scalar listZize(active_ib.size());
    label bodies_per_proc = ceil(listZize/Pstream::nProcs());
    InfoH << basic_Info << "Active IB listZize      : " << listZize<< endl;
    InfoH << basic_Info << "bodies_per_proc : " << bodies_per_proc<< endl;

    for(int assignProc = Pstream::myProcNo()*bodies_per_proc; assignProc < min((Pstream::myProcNo()+1)*bodies_per_proc,active_ib.size()); assignProc++)
    {
        const label body_id(active_ib[assignProc]);
        word path(curOutDir + "/body" + std::to_string(imm_bodies_[body_id].getBodyId()) +".info");
        OFstream ofStream(path);
        IOobject outClass
            (
                path,
                mesh_,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            );
        IOdictionary outDict(outClass);

        outDict.writeHeader(ofStream);
        imm_bodies_[body_id].record_body_info(outDict,curOutDir);
        outDict.writeData(ofStream);
    }

}
//---------------------------------------------------------------------------//
void IBMFoam::setSolverInfo()
{
    solver_info::setOnlyDEM(true);
}
//---------------------------------------------------------------------------//
