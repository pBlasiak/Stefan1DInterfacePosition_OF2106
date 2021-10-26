IOdictionary Stefan1DInterfacePositionDict
(
	IOobject
	(
		"Stefan1DInterfacePositionDict",
		runTime.constant(),
		mesh,
		IOobject::MUST_READ_IF_MODIFIED,
		IOobject::NO_WRITE
	)
);

// Phase change type
word phaseChangeType(Stefan1DInterfacePositionDict.lookup("phaseChangeType"));
// Automatic set of starting value for epsilon
bool epsilonStartingValue(Stefan1DInterfacePositionDict.lookupOrDefault("epsilonStartingValue", false));
// Wall temperature
dimensionedScalar Tw("Tw", dimTemperature, Stefan1DInterfacePositionDict);
