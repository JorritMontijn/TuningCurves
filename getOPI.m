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
	%2019-03-20 Created OSI function [by Jorrit Montijn]
	
	%check input
	if numel(vecTrialAngles) ~= size(matResp,2)
		error([mfilename ':NrOfTrialError'],'Number of trials in matrix does not match vector of angles');
	end
	
	% prep
	intN = size(matResp,1);
	vecOPI = nan(intN,1);
	[vecAngleIdx,vecUniqueAngles] = label2idx(vecTrialAngles);
	
	%run loop
	for intN=1:size(matResp,1)
		vecR = accumarray(vecAngleIdx',matResp(intN,:));
		vecOPI(intN) =  1 - circ_var(vecUniqueAngles(:),vecR);
	end
end

