function yOut = vonMisesDoubleFitBaseZero(xList,params)
	%vonMisesDoubleFit Creates output given x-data and parameters
	%	Syntax: yOut = vonMisesDoubleFit(xList,params)
	%	params must contain four values:
	%	- params(1): prefDir
	%	- params(2): kappa1
	%	- params(3): kappa2
	
	
	prefDir = params(1);
	kappa1 = params(2);
	kappa2 = params(3);
	
	if kappa1 < eps, kappa1 = eps;end
	if kappa2 < eps, kappa2 = eps;end
	
	[y1 x] = circ_vmpdf(xList, prefDir, kappa1);
	[y2 x] = circ_vmpdf(xList, (prefDir+pi), kappa2);
	
	yOut = y1+y2;
end

