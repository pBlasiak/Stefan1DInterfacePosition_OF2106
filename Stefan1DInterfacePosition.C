/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2008 OpenCFD Ltd.
     \\/     M anipulation  |
-------------------------------------------------------------------------------


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

// Interface position
dimensionedScalar X
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
    argList::addNote
    (
		"Calculates and writes the position of the interface for one dimensional"
		" Stefan problem."
		" It is assumed that the interface moves in the x direction."
    );

    argList::addOption("dict", "file", "Alternative setFieldsDict");

    #include "addRegionOption.H"

    timeSelector::addOptions();
    #include "setRootCase.H"
	#include "createTime.H"
    instantList timeDirs = timeSelector::select0(runTime, args);
	#include "createMesh.H"
	#include "createFields.H"
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

    Info<< "Number of iterations: " << iter    << endl;
	Info<< "epsilon = "             << epsilon << endl;

	// File for interface position
	OFstream IFfile("IFposition.txt");
	IFfile << "Time [s]\t" << "Numerical [m]\t" << "Analytical [m]\t" << "Error [%]" << endl;

	// File for velocity of the interface
	OFstream Ufile("IFvelocity.txt");
	Ufile << "Time [s]\t" << "Numerical [m]\t" << "Analytical [m]\t" << "Error [%]" << endl;

	// thermal diffusivity of vapor
    const dimensionedScalar D(k2/rho2/cp2);

	dimensionedScalar Uanal("Uanal", dimVelocity, 0.0);

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

		Info<< "Updating field U." << endl;
		volVectorField Uref
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

        Info<< endl;

		// Numerical value of interface position
        dimensionedScalar interfacePosition("interfacePosition", dimLength, 0.0);
	    interfacePosition = sum(mag(gradAlphal)*mesh.C().component(0))/max(sum(mag(gradAlphal)), dimensionedScalar("0gradAlpha", dimless/dimLength, SMALL));
		
		// Numerical value of interface velocity
		dimensionedScalar Unum("Unum", dimVelocity, 0.0);
		volVectorField cellsCenters(U.mesh().C());
		forAll(Uref, celli)
		{
			if
			(
				interfacePosition.value() >= cellsCenters[celli][0] 
			 && interfacePosition.value() <  cellsCenters[celli + 1][0]
			)
			{
				Unum.value() = Uref[celli][0];
			} 
		}

		if (phaseChangeType == "evaporation")
		{
			dimensionedScalar analyticalInterfacePosition= X(epsilon,k2, rho2, cp2, runTime.value());

			Info<<"Interface position for time " 
				<< runTime.timeName() 
				<< " is equal to: " 
				<< interfacePosition.value() 
				<< " [m]"
				<< " (error: "
				<< mag(interfacePosition.value() - analyticalInterfacePosition.value())/
						mag(analyticalInterfacePosition.value() + VSMALL)*100
				<< "%)"
				<< endl;

	    	Info<< "\nSaving the results to IFposition.txt\n" << endl;

			IFfile << runTime.timeName() 
	  		     << "\t" 
	  			 << interfacePosition.value() 
	  		     << "\t" 
	  			 << analyticalInterfacePosition.value() 
	  		     << "\t" 
				 << mag(interfacePosition.value() - analyticalInterfacePosition.value())/
						mag(analyticalInterfacePosition.value() + VSMALL)*100
	  			 << endl;

			if (runTime.value() == 0)
			{
				Ufile << runTime.timeName() 
	  			     << "\t" 
	  				 << Unum.value() 
	  			     << "\t" 
	  				 << 0 
	  			     << "\t" 
					 << mag(Unum.value() - Uanal.value())/
							mag(Uanal.value() + VSMALL)*100
	  				 << endl;
			}
			else
			{
				Uanal = (1.0 - rho2/rho1)*epsilon*Foam::sqrt(D/runTime);

	    		Info<< "\nSaving the results for interface velocity to IFvelocity.txt\n" << endl;

				Ufile << runTime.timeName() 
	  			     << "\t" 
	  				 << Unum.value() 
	  			     << "\t" 
	  				 << Uanal.value() 
	  			     << "\t" 
					 << mag(Unum.value() - Uanal.value())/
							mag(Uanal.value() + VSMALL)*100
	  				 << endl;
			}

			Info<< endl;
		}
		// ADD CALCULATIONS OF U FOR CONDENSATION!!!!!!!!!!!
		else
		{
			Info<<"Interface position for time " 
				<< runTime.timeName() 
				<< " is equal to: " 
				<< interfacePosition.value() 
				<< " [m]"
				<< " (error: "
				<< mag(interfacePosition.value() - X(epsilon, k1, rho1, cp1, runTime.value()).value())/
						mag(X(epsilon, k1, rho1, cp1, runTime.value()).value() + VSMALL)*100
				<< "%)"
				<< endl;

	    	Info<< "\nSaving the results to IFposition.txt\n" << endl;

			IFfile << runTime.timeName() 
	  		     << "\t" 
	  			 << interfacePosition.value() 
	  		     << "\t" 
	  			 << X(epsilon, k1, rho1, cp1, runTime.value()).value() 
	  		     << "\t" 
				 << mag(interfacePosition.value() - X(epsilon, k1, rho1, cp1, runTime.value()).value())/
						mag(X(epsilon, k1, rho1, cp1, runTime.value()).value() + VSMALL)*100
	  			 << endl;
		}

    }


	// Saves temperature distribution for the latest time
	word TfileName = "T_" + name(runTime.value()) + "s.txt";

	OFstream Tfile(TfileName);
	Tfile << "Time [s]\t"  << "x [m]\t" << "Numerical [m]\t" << "Analytical [m]\t" << "Error [%]" << endl;

	if (phaseChangeType == "evaporation")
	{
		volScalarField Tnum
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

		dimensionedScalar tau = runTime.value();
		dimensionedScalar avTau = k2/rho2/cp2*tau;
		dimensionedScalar analyticalInterfacePosition= X(epsilon,k2, rho2, cp2, runTime.value());
		
	   	Info<< "\nSaving the results for temperature at time " << name(runTime.value()) << " s to T_" << name(runTime.value()) << "s.txt\n" << endl;
		forAll(analyticalTemperature, celli)
		{
			x[celli] = analyticalTemperature.mesh().C()[celli].component(vector::X);
			if (x[celli] > analyticalInterfacePosition.value())
			{
				analyticalTemperature[celli] = TSat.value();
			}
			else
			{
				analyticalTemperature[celli] = Tw.value() + (TSat.value() - Tw.value())
					/Foam::erf(epsilon)*Foam::erf(x[celli]/2.0/Foam::sqrt(avTau.value())); 
			}
	
			Tfile << runTime.timeName() 
	  		     << "\t" 
	  			 << x[celli] 
	  		     << "\t" 
	  			 << Tnum[celli] 
	  		     << "\t" 
	  			 << analyticalTemperature[celli]
	  		     << "\t" 
				 << mag(Tnum[celli] - analyticalTemperature[celli])/
				 		analyticalTemperature[celli]*100
	  			 << endl;
		}
	}

	if (phaseChangeType == "condensation")
	{
		volScalarField Tnum
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

		dimensionedScalar tau = runTime.value();
		dimensionedScalar alTau = k1/rho1/cp1*tau;
		dimensionedScalar analyticalInterfacePosition= X(epsilon,k1, rho1, cp1, runTime.value());
		
	   	Info<< "\nSaving the results for temperature at time " << name(runTime.value()) << " s to T_" << name(runTime.value()) << "s.txt\n" << endl;
		forAll(analyticalTemperature, celli)
		{
			x[celli] = analyticalTemperature.mesh().C()[celli].component(vector::X);
			if (x[celli] > analyticalInterfacePosition.value())
			{
				analyticalTemperature[celli] = TSat.value();
			}
			else
			{
				analyticalTemperature[celli] = Tw.value() + (Tw.value() - TSat.value())
					/Foam::erf(epsilon)*Foam::erf(x[celli]/2.0/Foam::sqrt(alTau.value())); 
			}
	
			Tfile << runTime.timeName() 
	  		     << "\t" 
	  			 << x[celli] 
	  		     << "\t" 
	  			 << Tnum[celli] 
	  		     << "\t" 
	  			 << analyticalTemperature[celli]
	  		     << "\t" 
				 << mag(Tnum[celli] - analyticalTemperature[celli])/
				 		analyticalTemperature[celli]*100
	  			 << endl;
		}
	}

    Info<< "End" << endl;

    return 0;
}

// ************************************************************************* //
