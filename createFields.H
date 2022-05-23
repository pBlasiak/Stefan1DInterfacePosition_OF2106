Info<< "Reading field p_rgh\n" << endl;
volScalarField p_rgh
(
    IOobject
    (
        "p_rgh",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::NO_WRITE
    ),
    mesh
);
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
    p_rgh
);

volScalarField rho
(
    IOobject
    (
        "rho",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT
    ),
    mesh,
	dimensionedScalar("rho0", dimDensity, 0.0)
    //alpha1*rho1 + alpha2*rho2
);
rho.oldTime();

Info<< "Reading field T\n" << endl;
volScalarField T
(
    IOobject
    (
        "T",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::NO_WRITE
    ),
    mesh
);

Info << "Reading field U\n" << endl;
volVectorField U
(
    IOobject
    (
        "U",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::NO_WRITE
    ),
    mesh
);

Info<< "Creating field k\n" << endl;
volScalarField k
(
    IOobject
    (
        "k",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
	mesh,
    dimensionedScalar("k00", dimPower/dimLength/dimTemperature, 0)
);
// Need to store for turbulentHeatFluxTemperature BC
Info<< "Creating field cp\n" << endl;
volScalarField cp
(
    IOobject
    (
        "cp",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
	mesh,
    dimensionedScalar("cp00", dimensionSet(0,2,-2,-1,0,0,0), 0)
);
//cp = mixture->cp();

#   include "createPhi.H"

Info<< "Creating phaseChangeTwoPhaseMixture\n" << endl;
autoPtr<phaseChangeTwoPhaseMixture> mixture =
    phaseChangeTwoPhaseMixture::New(U, phi);

// liquid (1), vapor (2)
const dimensionedScalar& rho1 = mixture->rho1();
const dimensionedScalar& rho2 = mixture->rho2();
const dimensionedScalar& k1 = mixture->k1();
const dimensionedScalar& k2 = mixture->k2();
const dimensionedScalar& cp1 = mixture->cp1();
const dimensionedScalar& cp2 = mixture->cp2();
const dimensionedScalar& hEvap = mixture->hEvap();
const dimensionedScalar& TSat = mixture->TSatG();

Info << "Creating field gradAlphal\n" << endl;
volVectorField gradAlphal
(
    IOobject
    (
        "gradAlphal",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedVector("grad0", dimless/dimLength, vector(0,0,0))
);

volScalarField analyticalTemperature
(
    IOobject
    (
        "analyticalTemperature",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    TSat
);

volScalarField x
(
    IOobject
    (
        "analyticalTemperature",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("x0", dimLength, 0)
);
