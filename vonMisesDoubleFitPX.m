function yOut = vonMisesDoubleFitPX(params,xList)
	%vonMisesDoubleFit Creates output given x-data and parameters
	%	Syntax: yOut = vonMisesDoubleFit(xList,params)
	%	params must contain four values:
	%	- params(1): prefDir
	%	- params(2): kappa1
	%	- params(3): kappa2
	%	- params(4): baseline
	%	- params(5): gain
	
	prefDir = params(1);
	kappa1 = params(2);
	kappa2 = params(3);
	baseline = params(4);
	gain = params(5);
	
	if kappa1 < eps, kappa1 = eps;end
	if kappa2 < eps, kappa2 = eps;end
	
	[y1 x] = circ_vmpdf(xList, prefDir, kappa1);
	[y2 x] = circ_vmpdf(xList, (prefDir+pi), kappa2);
	
	yOut = gain*(y1+y2)+baseline;
end

