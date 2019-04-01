function [vecUniques,vecCounts] = getUniqueVals(vecValsIn)
	%getUniqueVals Retrieves a list of unique values from input vector
	%	Syntax: [vecUniques,vecCounts] = getUniqueVals(vecValsIn)
	%   Removes all duplicates from input list. If input is matrix, it
	%   automatically transforms the matrix to a 1D vector and performs the
	%   retrieval of unique values on the transformed vector
	%
	%Dependencies:
	% - none
	%
	%	Version history:
	%	1.0 - July 22 2013
	%	Created by Jorrit Montijn
	%	2.0 - Jan 10 2019
	%	Added counts [by JM]
	
	dblVal = -inf;
	vecUniques = [];
	vecCounts = [];
	dblMax = max(vecValsIn);
	vecValsIn = sort(vecValsIn(:),'ascend');
	while dblVal < dblMax
		[intInd] = find(vecValsIn > dblVal,1);
		dblVal = vecValsIn(intInd);
		vecUniques(end+1) = dblVal;
		vecCounts(end+1) = sum(vecValsIn==dblVal);
		vecValsIn = vecValsIn(vecValsIn > dblVal);
	end
end