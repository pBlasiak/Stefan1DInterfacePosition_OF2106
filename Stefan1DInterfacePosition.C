/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2008 OpenCFD Ltd.
     \\/     M anipulation  |
-------------------------------------------------------------------------------
2010-08-02 Eelco van Vliet: 1st public version of Stefan1DInterfacePosition:
  http://www.cfd-online.com/Forums/openfoam-solving/66705-wallheatflux-bc-not-constant-after-restart.html#post269812

2012-05-21 Eelco van Vliet:
  Quoting: http://www.cfd-online.com/Forums/openfoam-post-processing/101972-wallheatflux-utility-incompressible-case.html#post362191
  «modified the standard wallHeatflux utility which comes default with OF into
  a version for incompressible flows. Also removed a bug out of the code.»

2012-06-26 Eelco van Vliet:
  Quoting: http://www.cfd-online.com/Forums/openfoam-post-processing/101972-wallheatflux-utility-incompressible-case.html#post368330
  «p is now not required anymore.»

2014-06-22: Bruno Santos: Adapted to OpenFOAM 2.2.x.

-------------------------------------------------------------------------------
License
    This file is a derivative work of OpenFOAM.

    OpenFOAM is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 2 of the License, or (at your
    option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM; if not, write to the Free Software Foundation,
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

Application
    Stefan1DInterfacePosition

Description
    Calculates and writes the position of the interface for one dimensional
	Stefan problem. 
	It is assumed that the interface moves in the x direction.
\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
// modified from  wallHeatFlux
#include "singlePhaseTransportModel.H"
#include "phaseChangeTwoPhaseMixture.H"
#include "turbulentTransportModel.H"
#include "turbulenceModel.H"
#include "wallFvPatch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

scalar func(scalar x, const scalar constValue)
{
	return x*Foam::exp(pow(x,2))*Foam::erf(x) - constValue;
}

scalar derivErf(const scalar x)
{
	return 2*Foam::exp(-pow(x,2))/Foam::sqrt(Foam::constant::mathematical::pi);
}

scalar derivFunc(scalar x)
{
	return 
		Foam::exp(pow(x,2))*Foam::erf(x) 
	  + x*(2*x*Foam::exp(pow(x,2))*Foam::erf(x) + Foam::exp(pow(x,2))*derivErf(x) ); 
}

dimensionedScalar deltaInter
       (
	       const scalar eps, 
		   const dimensionedScalar thermCond, 
		   const dimensionedScalar rho, 
		   const dimensionedScalar cp, 
		   const scalar time
	   )
{
	return 2*eps*Foam::sqrt(thermCond*time/rho/cp);
}

int main(int argc, char *argv[])
{

    timeSelector::addOptions();
    #include "setRootCase.H"
#   include "createTime.H"
    instantList timeDirs = timeSelector::select0(runTime, args);
#   include "createMesh.H"
#       include "createFields.H"
	#include "readTransportProperties.H"

	const label maxIterNr = 20000;
	label iter = 0;
	const scalar tol = 1e-12;
	dimensionedScalar LHS;
	scalar epsilon = 0.1; // guess of starting value
	scalar epsilonPrev = 0;

	if (phaseChangeType == "evaporation")
	{
		LHS = cp2*(Tw - TSat)/hEvap/Foam::sqrt(Foam::constant::mathematical::pi);
	}
	else if (phaseChangeType == "condensation")
	{
		LHS = cp1*(TSat - Tw)/hEvap/Foam::sqrt(Foam::constant::mathematical::pi);
	}
	else
	{
		FatalErrorIn("Stefan1DInterfacePosition") 
			<< "phaseChangeType can be \"evaporation\" or \"condensation\" only"
			<< nl
			<< exit(FatalError);
	}

	if (epsilonStartingValue)  
	{
	    Info << "Type the starting value for epsilon:" << endl;	
	    std::cin >> epsilon;
	}
	

	do 
	{
		epsilonPrev = epsilon;		
		epsilon -= func(epsilon, LHS.value())/derivFunc(epsilon);	
		iter++;
	} while ( mag(epsilon - epsilonPrev) > tol && iter < maxIterNr );

    Info<< "Liczba iteracji: " << iter << endl;
	Info<< "epsilon = " << epsilon << endl;

	OFstream IFfile("IFposition.txt");
	IFfile << "Time [s]\t" << "Numerical [m]\t" << "Analytical [m]\t" << "Error [%]" << endl;

    forAll(timeDirs, timeI)
    {
        runTime.setTime(timeDirs[timeI], timeI);
        Info<< "Time = " << runTime.timeName() << endl;
        mesh.readUpdate();

		Info<< "Reading field alpha." << mixture->phase1Name() << endl;
    	volScalarField alphal
    	(
    	    IOobject
    	    (
				"alpha." + mixture->phase1Name(),
    	        runTime.timeName(),
    	        mesh,
    	        IOobject::MUST_READ,
    	        IOobject::NO_WRITE
    	    ),
    	    mesh
    	);

        gradAlphal=fvc::grad(alphal);

        Info<< endl;

        dimensionedScalar interfacePosition("interfacePosition", dimLength, 0.0);
   
	    interfacePosition = sum(mag(gradAlphal)*mesh.C().component(0))/sum(mag(gradAlphal));


		if (phaseChangeType == "evaporation")
		{
			Info<<"Interface position for time " 
				<< runTime.timeName() 
				<< " is equal to: " 
				<< interfacePosition.value() 
				<< " [m]"
				<< " (error: "
				<< mag(interfacePosition.value() - deltaInter(epsilon,k2, rho2, cp2, runTime.value()).value())/
						mag(deltaInter(epsilon, k2, rho2, cp2, runTime.value()).value() + VSMALL)*100
				<< "%)"
				<< endl;

	    	Info<< "\nSaving the results to IFposition.txt\n" << endl;

			IFfile << runTime.timeName() 
	  		     << "\t" 
	  			 << interfacePosition.value() 
	  		     << "\t" 
	  			 << deltaInter(epsilon, k2, rho2, cp2, runTime.value()).value() 
	  		     << "\t" 
				 << mag(interfacePosition.value() - deltaInter(epsilon,k2, rho2, cp2, runTime.value()).value())/
						mag(deltaInter(epsilon, k2, rho2, cp2, runTime.value()).value() + VSMALL)*100
	  			 << endl;
		}
		else
		{
			Info<<"Interface position for time " 
				<< runTime.timeName() 
				<< " is equal to: " 
				<< interfacePosition.value() 
				<< " [m]"
				<< " (error: "
				<< mag(interfacePosition.value() - deltaInter(epsilon, k1, rho1, cp1, runTime.value()).value())/
						mag(deltaInter(epsilon, k1, rho1, cp1, runTime.value()).value() + VSMALL)*100
				<< "%)"
				<< endl;

	    	Info<< "\nSaving the results to IFposition.txt\n" << endl;

			IFfile << runTime.timeName() 
	  		     << "\t" 
	  			 << interfacePosition.value() 
	  		     << "\t" 
	  			 << deltaInter(epsilon, k1, rho1, cp1, runTime.value()).value() 
	  		     << "\t" 
				 << mag(interfacePosition.value() - deltaInter(epsilon, k1, rho1, cp1, runTime.value()).value())/
						mag(deltaInter(epsilon, k1, rho1, cp1, runTime.value()).value() + VSMALL)*100
	  			 << endl;
		}
    }

    Info<< "End" << endl;

    return 0;
}

// ************************************************************************* //
