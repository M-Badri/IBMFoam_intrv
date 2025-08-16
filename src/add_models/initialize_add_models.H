    // Note: if you want to add a new add_model, do it here

    scalar thrSurf(readScalar(IBMFoamDict_.lookup("surface_threshold")));

    forAll (body_names_,nameI)
    {
        word body_name(body_names_[nameI]);
        InfoH << add_model_Info << "-- add_modelMessage-- "
            << "Creating addition model for body: " << body_name << endl;

       
        dictionary add_modelDict(IBMFoamDict_.subDict(body_name).subDict("bodyAddition"));
        word howToAdd(word(add_modelDict.lookup("add_model")));
        word b_geom;
        std::unique_ptr<geom_model> b_geomModel;
      
        if (IBMFoamDict_.subDict(body_name).found("b_geom"))
        {
            word input = IBMFoamDict_.subDict(body_name).lookup("b_geom");
            b_geom = input;
            InfoH << iB_Info << "Found b_geom for " << body_name
                << ", the body is: " << b_geom << endl;
        }
        else
        {
            b_geom = "convex";
            InfoH << iB_Info << "Did not find b_geom for "
                << body_name << ", using b_geom: " << b_geom << endl;
        }

        if(b_geom == "convex")
        {
            word stlPath("constant/triSurface/"+body_name+".stl");
            b_geomModel = std::unique_ptr<convex_body>(
                new convex_body(mesh_,stlPath,thrSurf));
        }
        else if(b_geom == "nonConvex")
        {
            word stlPath("constant/triSurface/"+body_name+".stl");
            b_geomModel = std::unique_ptr<nonConvex_body>(
                new nonConvex_body(mesh_,stlPath,thrSurf));
        }
        else if(b_geom == "sphere")
        {
            vector startPosition = IBMFoamDict_.subDict(body_name).subDict("sphere").lookup("startPosition");
            scalar radius = readScalar(IBMFoamDict_.subDict(body_name).subDict("sphere").lookup("radius"));

            b_geomModel = std::unique_ptr<sphere_body>(
                new sphere_body(mesh_,startPosition,radius,thrSurf));
        }
        else
        {
            word stlPath("constant/triSurface/"+body_name+".stl");
            InfoH << iB_Info << "b_geom: " << b_geom
                << " not supported, using b_geom nonConvex" << endl;
            b_geom = "nonConvex";
            b_geomModel = std::unique_ptr<nonConvex_body>(
                new nonConvex_body(mesh_,stlPath,thrSurf));
        }

        InfoH << add_model_Info << "-- add_modelMessage-- "
            << "trying to set add_model: " << howToAdd << endl;

        if (howToAdd == "once")
        {
            add_models_.set
            (
                nameI,
                new add_modelOnce
                (
                    add_modelDict,
                    mesh_,
                    startTime0,
                    std::move(b_geomModel),
                    cellPoints_
                )
            );
        }
        else if (howToAdd == "onceFromFile")
        {
            add_models_.set
            (
                nameI,
                new add_modelOnceFromFile
                (
                    add_modelDict,
                    mesh_,
                    startTime0,
                    std::move(b_geomModel),
                    cellPoints_,
                    b_geom,
                    thrSurf
                )
            );
        }
        else if (howToAdd == "onceScatter")
        {
            add_models_.set
            (
                nameI,
                new add_modelOnceScatter
                (
                    add_modelDict,
                    mesh_,
                    startTime0,
                    std::move(b_geomModel),
                    cellPoints_
                )
            );
        }
        else if (howToAdd == "repeatSamePosition")
        {
            add_models_.set
            (
                nameI,
                new add_modelRepeatSamePosition
                (
                    add_modelDict,
                    mesh_,
                    std::move(b_geomModel),
                    cellPoints_
                )
            );
        }
        else if (howToAdd == "repeatRandomPosition")
        {
            add_models_.set
            (
                nameI,
                new add_modelRepeatRandomPosition
                (
                    add_modelDict,
                    mesh_,
                    std::move(b_geomModel),
                    cellPoints_
                )
            );
        }
        else if (howToAdd == "distribution")
        {
            add_models_.set
            (
                nameI,
                new add_modelDistribution
                (
                    add_modelDict,
                    mesh_,
                    std::move(b_geomModel),
                    cellPoints_
                )
            );
        }
        else
        {
            InfoH << add_model_Info << "-- add_modelMessage-- "
                << "No addition model specified, defaulting to add_modelOnce"
                << endl;

            add_models_.set
            (
                nameI,
                new add_modelOnce
                (
                    add_modelDict,
                    mesh_,
                    startTime0,
                    std::move(b_geomModel),
                    cellPoints_
                )
            );
        }
    }
