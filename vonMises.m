function yOut = vonMises(xList,params)
	%vonMises Creates output given x-data and parameters
	%	Syntax: yOut = vonMises(xList,params)
	%	params must contain four values:
	%	- params(1): prefDir
	%	- params(2): kappa
	%	- params(3): baseline
	%	- params(4): scale factor
	
	prefDir = params(1);
	kappa = params(2);
	baseline = 0;%params(3);
	dblScaleFactor = 1;%params(4);
	
	if kappa < eps, kappa = eps;end
	
	[y x] = circ_vmpdf(xList, prefDir, kappa);

	yOut = imnorm(y');
end