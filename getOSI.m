function [vecOSI,vecAngleResp,intIdx,intOrthIdx] = getOSI(matResp,vecTrialAngles)
	%getOSI Calculates the orientation selectivity index
	%	 vecOSI = getOSI(matResp,vecTrialAngles)
	%
	%OSI = (Rpref_ori - Rorth)/Rpref_ori
	%
	%Inputs:
	% - matResp; [Neuron x Trial] response vector
	% - vecAngles; [1 x Trial] stimulus orientation vector  [in radians]
	%
	%Version History:
	%2019-03-20 Created OSI function [by Jorrit Montijn]
	
	if numel(vecTrialAngles) ~= size(matResp,2)
		error([mfilename ':NrOfTrialError'],'Number of trials in matrix does not match vector of angles');
	end
	
	% prep
	matResp = bsxfun(@minus,matResp,min(matResp,[],2));
	intN = size(matResp,1);
	vecOSI = nan(intN,1);
	[vecAngleIdx,vecUniqueAngles] = label2idx(vecTrialAngles);
	intStimNum = numel(vecUniqueAngles);
	
	%run loop
	for intN=1:size(matResp,1)
		vecTrialResp = matResp(intN,:);
		vecTrialAngleIdx = vecAngleIdx;
		vecTrialAngleIdx(isnan(vecTrialResp)) = [];
		vecTrialResp(isnan(vecTrialResp)) = [];
		vecAngleResp = nan(1,intStimNum);
		for intStimIdx=1:intStimNum
			vecAngleResp(intStimIdx) = nanmean(vecTrialResp(vecTrialAngleIdx==intStimIdx));
		end
		
		% get OI
		[dblPref,intIdx] = max(vecAngleResp);
		dblPrefAngle = vecUniqueAngles(intIdx);
		[dblOrthAngle,intOrthIdx] = min(abs(circ_dist(vecUniqueAngles,dblPrefAngle+deg2rad(90))));
		if intOrthIdx == 0,intOrthIdx=numel(vecUniqueAngles);end
		vecOSI(intN) = (dblPref - vecAngleResp(intOrthIdx)) ./ dblPref;
	end
	
end

