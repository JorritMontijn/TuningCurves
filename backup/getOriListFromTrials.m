function vecOriOut = getOriListFromTrials(vecOriIn)
	intOri = -inf;
	vecOriOut = [];
	intMax = max(vecOriIn);
	vecOriIn = sort(vecOriIn,'ascend');
	while intOri < intMax
		[intInd] = find(vecOriIn > intOri,1);
		intOri = vecOriIn(intInd);
		vecOriIn = vecOriIn(vecOriIn > intOri);
		vecOriOut(end+1) = intOri;
	end
end