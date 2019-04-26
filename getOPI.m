function vecOPI = getOPI(matResp,vecTrialAngles)
	%getOPI Calculates the orientation precision index (Ringach)
	%	 vecOPI = getOPI(matResp,vecTrialAngles)
	%
	%OPI = 1 - circ_var
	%
	%Inputs:
	% - matResp; [Neuron x Trial] response vector
	% - vecAngles; [1 x Trial] stimulus orientation vector  [in radians]
	%
	%Version History:
	%2019-03-20 Created OPI function [by Jorrit Montijn]
	
	%check input
	if numel(vecTrialAngles) ~= size(matResp,2)
		error([mfilename ':NrOfTrialError'],'Number of trials in matrix does not match vector of angles');
	end
	
	% prep
	matResp = matResp - min(matResp,[],2);
	intN = size(matResp,1);
	vecOPI = nan(intN,1);
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
		%vecAngleResp = vecAngleResp - min(vecAngleResp);
		vecOPI(intN) =  1 - circ_var(vecUniqueAngles(:),vecAngleResp(:));
	end
end

