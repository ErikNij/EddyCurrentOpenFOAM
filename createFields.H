// Block vector field for magnetic vector potential real part (first entry) and imaginary part (second entry)
//Not needed
//\Not needed
    Info<< "Reading dictionary: physicalProperties\n" << endl;

// read .../constant/physicalerProperties
    IOdictionary physicalProperties
    (
        IOobject
        (
            "physicalProperties",
            runTime.constant(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    );

    dimensionedScalar muMag
    (
        physicalProperties.lookup("muMag")
    );


    dimensionedScalar w
    (
        physicalProperties.lookup("w")
    );

    volScalarField sigma
    (
        IOobject
        (
            "sigma",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );

    Info<< "Reading intial conditions of field variables from .../0 directory\n" << endl;

    volVectorField JCoil
    (
        IOobject
        (
            "JCoil",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );

    volVectorField HCoil
    (
        IOobject
        (
            "HCoil",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedVector("HCoil", dimensionSet(0,0,0,0,0), Foam::Vector<double>(0.0, 0.0, 0.0)) // All zeros
    );

    //Debug suite: Here are a bunch of outputs that are useful for debugging:


    volScalarField EdgeCell
        (
        IOobject
            (
                "EdgeCell",  // Field name
                runTime.timeName(),
                mesh,            
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            mesh,
            dimensionedScalar("EdgeCell", dimensionSet(0,0,0,0,0), 0) // All zeros
    );

    volScalarField divJI
        (
        IOobject
            (
                "divJI",  // Field name
                runTime.timeName(),
                mesh,            
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            mesh,
            dimensionedScalar("divJI", dimensionSet(0,0,0,0,0), 0) // All zeros
    );


    volVectorField RestrictVec
        (
        IOobject
            (
                "RestrictVec",  // Field name
                runTime.timeName(),
                mesh,            
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh,
            dimensionedVector("RestrictVec", dimensionSet(0,0,0,0,0), Foam::Vector<double>(0.0, 0.0, 0.0)) // All zeros
    );

    //volVectorField HCoil
    //    (
    //    IOobject
    //        (
    //            "HCoil",  // Field name
    //            runTime.timeName(),
    //            mesh,            
    //            IOobject::MUST_READ,
    //            IOobject::AUTO_WRITE
    //        ),
    //        mesh,
    //        dimensionedVector("HCoil", dimensionSet(0,0,0,0,0), Foam::Vector<double>(0.0, 0.0, 0.0)) // All zeros
    //);


    volVectorField TR
        (
        IOobject
            (
                "TR",  // Field name
                runTime.timeName(),
                mesh,            
                IOobject::MUST_READ,
                IOobject::AUTO_WRITE
            ),
            mesh
    );

    volVectorField TI
        (
        IOobject
            (
                "TI",  // Field name
                runTime.timeName(),
                mesh,            
                IOobject::MUST_READ,
                IOobject::AUTO_WRITE
            ),
            mesh
    );

    volVectorField TINew
        (
        IOobject
            (
                "TINew",  // Field name
                runTime.timeName(),
                mesh,            
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh,
            dimensionedVector("TINew", dimensionSet(0,0,0,0,0), Foam::Vector<double>(0.0, 0.0, 0.0)) // All zeros
    );

    volVectorField TRNew
        (
        IOobject
            (
                "TRNew",  // Field name
                runTime.timeName(),
                mesh,            
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh,
            dimensionedVector("TRNew", dimensionSet(0,0,0,0,0), Foam::Vector<double>(0.0, 0.0, 0.0)) // All zeros
    );


    //volScalarField FitnessOld
    //    (
    //    IOobject
    //        (
    //            "FitnessOld",  // Field name
    //            runTime.timeName(),
    //            mesh,            
    //            IOobject::NO_READ,
    //            IOobject::AUTO_WRITE
    //        ),
    //        mesh,
    //        dimensionedScalar("FitnessOld", dimensionSet(0,0,0,0,0), 0) // All zeros
    //);
//
    //volScalarField FitnessNew
    //    (
    //    IOobject
    //        (
    //            "FitnessNew",  // Field name
    //            runTime.timeName(),
    //            mesh,            
    //            IOobject::NO_READ,
    //            IOobject::AUTO_WRITE
    //        ),
    //        mesh,
    //        dimensionedScalar("FitnessNew", dimensionSet(0,0,0,0,0), 0) // All zeros
    //);
//
    
    

    //volScalarField Tmag
    //(
    //    IOobject
    //    (
    //        "Tmag",
    //        runTime.timeName(),
    //        mesh,
    //        IOobject::NO_READ,
    //        IOobject::AUTO_WRITE
    //    ),
    //    sqrt((TR & TR) + (TI & TI)) // T = sqrt(TR*TR + TI*TI) 
    //);
//


    volVectorField HR
        (
        IOobject
            (
                "HR",  // Field name
                runTime.timeName(),
                mesh,            
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            mesh,
            dimensionedVector("HR", dimensionSet(0,0,0,0,0), Foam::Vector<double>(0.0, 0.0, 0.0)) // All zeros
    );

    

    volVectorField HIndR
        (
        IOobject
            (
                "HIndR",  // Field name
                runTime.timeName(),
                mesh,            
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            mesh,
            dimensionedVector("HIndR", dimensionSet(0,0,0,0,0), Foam::Vector<double>(0.0, 0.0, 0.0)) // All zeros
    );

    volVectorField HIndI
        (
        IOobject
            (
                "HIndI",  // Field name
                runTime.timeName(),
                mesh,            
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            mesh,
            dimensionedVector("HIndI", dimensionSet(0,0,0,0,0), Foam::Vector<double>(0.0, 0.0, 0.0)) // All zeros
    );

    volVectorField HI
        (
        IOobject
            (
                "HI",  // Field name
                runTime.timeName(),
                mesh,            
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            mesh,
            dimensionedVector("HI", dimensionSet(0,0,0,0,0), Foam::Vector<double>(0.0, 0.0, 0.0)) // All zeros
    );

//
    volScalarField Hmag
    (
        IOobject
        (
            "Hmag",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        sqrt((HR & HR) + (HI & HI))     
    );
//
    volScalarField Bmag
    (
        IOobject
        (
            "Bmag",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        Hmag * muMag
    );
//
    volVectorField JR
        (
            IOobject
            (
            "JR",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        fvc::curl(TR)
    );
//
    volVectorField JI
        (
            IOobject
            (
            "JI",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        fvc::curl(TI)
    );
//
    //volScalarField Jind
    //(
    //    IOobject
    //    (
	//		"Jind",
	//		runTime.timeName(),
	//		mesh,
	//		IOobject::NO_READ,
	//		IOobject::AUTO_WRITE
    //    ),
    //    sqrt((JR & JR) + (JI & JI))   
    //);
//
	//volScalarField qJ
	//(
	//	IOobject
	//	(
	//		"qJ",
	//		runTime.timeName(),
	//		mesh,
	//		IOobject::NO_READ,
	//		IOobject::AUTO_WRITE
	//	),
	//	sigma/2.0 * sqr(w) * sqr(Tmag)  
	//);
//
    //volVectorField fL
    //    (
    //        IOobject
    //        (
    //        "fL",
    //        runTime.timeName(),
    //        mesh,
    //        IOobject::NO_READ,
    //        IOobject::AUTO_WRITE
    //    ),
    //    0.5 * ((JR ^ (HR * muMag)) - (JI ^ (HI * muMag)))
    //);
//
//