if (mesh_.nGeometricD() < 3)
    {
        const vector validDirs = (mesh_.geometricD() + Vector<label>::one)/2;
        Axis_ -= validDirs;
    }
    Axis_ /= mag(Axis_);

    if (immersedDict_.found("staticBody"))
    {
        bodyOperation_ = 0;
        InfoH << iB_Info << body_name_ << " is a static body." << endl;
    }
    else if (immersedDict_.found("prescribedTransBody"))
    {
        bodyOperation_ = 1;

        Vel_   = immersedDict_.subDict("prescribedTransBody").lookup("velocity");

        InfoH << iB_Info << body_name_ << " is a freely rotating body with prescribed linear velocity." << endl;
    }
    else if (immersedDict_.found("prescribedRotBody"))
    {
        bodyOperation_ = 2;

        Axis_ = immersedDict_.subDict("prescribedRotBody").lookup("axis");
        omega_  = readScalar(
                    immersedDict_.subDict("prescribedRotBody").lookup("omega")
                );

        InfoH << iB_Info << body_name_ << " is a freely moving body with prescribed rotation." << endl;
    }
    else if (immersedDict_.found("prescribedTransRotBody"))
    {
        bodyOperation_ = 3;

        Vel_   = immersedDict_.subDict("prescribedTransRotBody").lookup("velocity");
        Axis_ = immersedDict_.subDict("prescribedTransRotBody").lookup("axis");
        omega_  = readScalar(
                    immersedDict_.subDict("prescribedTransRotBody").lookup("omega")
                );

        InfoH << iB_Info << body_name_ << " has prescribed both movement and rotation." << endl;
    }
    else if (immersedDict_.found("prescribedTransFixedAxisRotBody"))
    {
        bodyOperation_ = 4;

        Vel_   = immersedDict_.subDict("prescribedTransFixedAxisRotBody").lookup("velocity");
        Axis_ = immersedDict_.subDict("prescribedTransFixedAxisRotBody").lookup("axis");

        InfoH << iB_Info << body_name_ << " has prescribed movement and axis of rotation." << endl;
    }
    else if (immersedDict_.found("fullyCoupledBody"))
    {
        bodyOperation_ = 5;

        if (immersedDict_.isDict("fullyCoupledBody"))
        {
            if(immersedDict_.subDict("fullyCoupledBody").found("velocity"))
            {
                Vel_ = immersedDict_.subDict("fullyCoupledBody").lookup("velocity");
            }
            if(immersedDict_.subDict("fullyCoupledBody").found("axis"))
            {            
                Axis_ = immersedDict_.subDict("fullyCoupledBody").lookup("axis");
            }
            if(immersedDict_.subDict("fullyCoupledBody").found("omega"))
            {
                omega_  = readScalar(
                    immersedDict_.subDict("fullyCoupledBody").lookup("omega")
                );
            }
        }

        InfoH << iB_Info << body_name_ << " is fully coupled with fluid phase." << endl;
    }
    else
    {
        InfoH << iB_Info << "No body operation was found for " << body_name_ << endl
          << "Assuming static body.";
    }

    if (immersedDict_.found("update_torq"))
    {
        update_torq_ = readBool(immersedDict_.lookup("update_torq"));
        InfoH << iB_Info << "Found update_torq: " << update_torq_ << endl;
    }
    else
    {
        InfoH << iB_Info << "Did not find update_torq, using update_torq: " << update_torq_ << endl;
    }

    if (immersedDict_.found("startSynced"))
    {
        startSynced_ = readBool(immersedDict_.lookup("startSynced"));
        if (startSynced_)
        {
            InfoH << iB_Info << "Will try to sync the body with the flow upon creation" << endl;
        }
        else
        {
            InfoH << iB_Info << "The body will be created as static" << endl;
        }
    }
    else
    {
        InfoH << iB_Info << "startSynced was not specified, using startSynced: " << startSynced_ << endl;
    }

    if (immersedDict_.found("sdBasedLambda"))
    {
        sdBasedLambda_ = readBool(immersedDict_.lookup("sdBasedLambda"));
    }
    if (immersedDict_.found("interface_span"))
    {
        intSpan_ = readScalar(immersedDict_.lookup("interface_span"));
    }
    if (immersedDict_.found("refineBuffers"))
    {
        refineBuffers_ = readLabel(immersedDict_.lookup("refineBuffers"));
    }
    if (immersedDict_.found("t_to_set_static"))
    {
        t_to_set_static_ = readLabel(immersedDict_.lookup("t_to_set_static"));
    }

    geom_model_->setIntSpan(intSpan_);
    geom_model_->setSdBasedLambda(sdBasedLambda_);
    geom_model_->setRhoS(immersedDict_.lookup("rho"));

    string material = word(immersedDict_.lookup("material"));

    if(ibIntp.valid())
    {
        intpInfo_.set(ibIntp().createIntpInfo(
            mesh_,
            geom_model_
        ));
    }

    ibContactClass_.set(new ibContactClass(
        geom_model_,
        material
    ));

    ibCntVars_.set(new ibContactVars
    (
        body_id_,
        Vel_,
        omega_,
        Axis_,
        geom_model_->getM0(),
        geom_model_->getM(),
        geom_model_->getRhoS()
    ));

    wallCntInfo_.set(new wallContactInfo
    (
        ibContactClass_(),
        ibCntVars_()
    ));
