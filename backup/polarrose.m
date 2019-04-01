function h = polarrose(theta,rho)
	%POLARROSE   Angle histogram plot.
	%   POLARROSE(THETA,RHO) plots the angle histogram for the angles in
	%   THETA with the magnitude being RHO
	%   H = ROSE(...) returns the figure handle.
	%
	%   See also ROSE, HIST, POLAR, COMPASS.
	
	%   Clay M. Thompson 7-9-91
	%
	%	Jorrit S. Montijn 4-9-12
	%	Transformed function to work with theta, rho input syntax
	
	% Determine bin edges and get histogram
	edges = sort(rem([(theta(2:end)+theta(1:end-1))/2 (theta(end)+theta(1)+2*pi)/2],2*pi));
	edges = [edges edges(1)+2*pi];
	
	nn = circshift(rho,[0 -1]);
	
	% Form radius values for histogram triangle
	if min(size(nn))==1, % Vector
		nn = nn(:);
	end
	[m,n] = size(nn);
	mm = 4*m;
	r = zeros(mm,n);
	r(2:4:mm,:) = nn;
	r(3:4:mm,:) = nn;
	
	% Form theta values for histogram triangle from triangle centers (xx)
	zz = edges;
	
	t = zeros(mm,1);
	t(2:4:mm) = zz(1:m);
	t(3:4:mm) = zz(2:m+1);
	
	h = polar(t,r);
end