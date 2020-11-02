function yOut = vonMisesDoubleFitPX(params,xList)
	%vonMisesDoubleFit Creates output given x-data and parameters
	%	Syntax: yOut = vonMisesDoubleFit(xList,params)
	%	params must contain four values:
	%	- params(1): prefDir
	%	- params(2): kappa
	%	- params(3): direction index
	%	- params(4): baseline
	%	- params(5): gain
	
	prefDir = params(1);
	kappa = params(2);
	dirIdx = params(3);
	baseline = params(4);
	gain = params(5);
	
	if kappa < eps, kappa = eps;end
	
	[y1 x] = circ_vmpdf(xList, prefDir, kappa);
	[y2 x] = circ_vmpdf(xList, (prefDir+pi), kappa);
	
	yOut = reshape(gain*(y1+dirIdx*y2)+baseline,size(xList));
end

