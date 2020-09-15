function yOut = vonMisesSingleFitPX(params,xList)
	%vonMisesSingleFitPX Creates output given parameters and x-data
	%	Syntax: yOut = vonMisesSingleFitPX(params,xList)
	%	params must contain four values:
	%	- params(1): prefDir
	%	- params(2): kappa
	%	- params(3): baseline
	%	- params(4): gain
	
	prefDir = params(1);
	kappa = params(2);
	baseline = params(3);
	gain = params(4);
	
	if kappa < eps, kappa = eps;end
	
	[y x] = circ_vmpdf(xList, prefDir, kappa);
	
	yOut = gain*reshape(y,size(xList))+baseline;
end
