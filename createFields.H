Info<< "Reading field p_rgh\n" << endl;
volScalarField p_rgh
(
    IOobject
    (
        "p_rgh",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

Info<< "Reading field U\n" << endl;
volVectorField U
(
    IOobject
    (
        "U",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

#include "createPhi.H"


Info<< "Reading incompressibleTwoPhaseInteractingMixture\n" << endl;
incompressibleTwoPhaseInteractingMixture mixture(U, phi);

volScalarField& alpha1(mixture.alpha1());


// Mixture density
volScalarField rho
(
    IOobject
    (
        "rho",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mixture.rho()
);

// Mixture viscosity
volScalarField muMix
(
    IOobject
    (
        "muMix",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mixture.mu()
);

volScalarField strRate
(
    IOobject
    (
        "strRate",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    Foam::sqrt(2.0)*mag(symm(fvc::grad(U)))
);

// Mass flux
surfaceScalarField rhoPhi
(
    IOobject
    (
        "rhoPhi",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    fvc::interpolate(rho)*phi
);


// Relative Velocity
autoPtr<relativeVelocityModel> UdmModelPtr
(
    relativeVelocityModel::New
    (
        mixture,
        mixture
    )
);


// Construct compressible turbulence model
autoPtr
<
    CompressibleTurbulenceModel<incompressibleTwoPhaseInteractingMixture>
> turbulence
(
    CompressibleTurbulenceModel<incompressibleTwoPhaseInteractingMixture>
    ::New(rho, U, rhoPhi, mixture)
);


#include "readGravitationalAcceleration.H"
#include "readhRef.H"
#include "gh.H"


volScalarField p
(
    IOobject
    (
        "p",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    p_rgh + rho*gh
);

label pRefCell = 0;
scalar pRefValue = 0.0;
setRefCell
(
    p,
    p_rgh,
    pimple.dict(),
    pRefCell,
    pRefValue
);

if (p_rgh.needReference())
{
    p += dimensionedScalar
    (
        "p",
        p.dimensions(),
        pRefValue - getRefCellValue(p, pRefCell)
    );
    p_rgh = p - rho*gh;
}

mesh.setFluxRequired(p_rgh.name());
mesh.setFluxRequired(alpha1.name());

// MULES Correction
tmp<surfaceScalarField> talphaPhiCorr0;

#include "createMRF.H"


//-Adding properties for the tracer equations to be solved

Info<< "Reading field C\n" << endl;

    volScalarField C
    (
        IOobject
        (
            "C",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );

    Info<< "Reading diffusivity\n" << endl;

    volScalarField DT
    (
        IOobject
        (
            "DT",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );


    Info<< "Reading transportProperties\n" << endl;

    IOdictionary transportProperties
    (
        IOobject
        (
            "transportProperties",
            runTime.constant(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    );

    Info<< "Reading Schmidt Number\n" << endl;

    dimensionedScalar Sc
    (
        transportProperties.lookup("Sc")
    );


    Info<< "Reading Molecular Diffusivity\n" << endl;

    dimensionedScalar Dm
    (
        transportProperties.lookup("Dm")
    );


//------------------------Compression Terms
    volScalarField dComp
    (
        IOobject
        (
            "dComp",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
	dimensionedScalar("dComp", dimensionSet(0,2,-1,0,0,0,0), 0.0)		//-initial value
    );

    Info<< "Reading betaComp\n" << endl;	//-unit must be kg/m3

    dimensionedScalar betaComp
    (
        transportProperties.lookup("betaComp")
    );

   Info<< "Reading lamdaComp\n" << endl;	//-unit must be kg/ms2

    dimensionedScalar lambdaComp
    (
        transportProperties.lookup("lambdaComp")
    );

    Info<< "Reading Critical volume-fraction\n" << endl;	//-unit must be unitless

    dimensionedScalar alphaComp
    (
        transportProperties.lookup("alphaComp")
    );
    //----magnitude of g for calculating the compression term
    dimensionedScalar gComp
    (
        transportProperties.lookup("gComp")
    );

   //----Dummy term so that I keep OF happy with the dimensions
    dimensionedScalar dummyComp
    (
        transportProperties.lookup("dummyComp")
    );



