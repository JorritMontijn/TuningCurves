function [matRespNSR,vecStimTypes,vecUnique] = getStimulusResponses(matResp,vecStimTypeList)
	%getStimulusResponses Transforms [Neuron x Trial] to [Neuron x Stim x Rep]
	%Syntax:
	%[matRespNSR,vecStimTypes,vecUnique] = getStimulusResponses(matResp,vecStimTypeList)
	%
	%Inputs:
	% - matResp;		response matrix [Neuron x Trial] 
	% - vecStimTypeList;stimulus vector (e.g., angle or ID) [1 x Trial] 
	%
	%Outputs:
	% - matRespNSR;		response matrix [Neuron x Stim x Rep]
	% - vecStimTypes;	indexed values corresponding to vecStimTypeList
	% - vecUnique;		unique stimulus values from vecStimTypeList that
	%					correspond to the indexed values in vecStimTypes
	%
	%	By Jorrit Montijn (Alex Pouget lab), 22-02-18 (dd-mm-yy; Universite de Geneve)
	%%
	%neuron x trial
	intNeurons = size(matResp,1);
	intTrials = size(matResp,2);
	[vecStimTypes,vecUnique,vecRepetitions] = val2idx(vecStimTypeList);
	vecStimTypes = vecStimTypes(:)';
	intStimTypes = numel(unique(vecStimTypes));
	intRepetitions = max(vecRepetitions);
	if ~all(vecRepetitions==intRepetitions)
		warning([mfilename ':InconsistentNumberOfRepetitions'],'Number of repetitions is inconsistent');
	end
	
	matRespNSR = nan(intNeurons,intStimTypes,intRepetitions);
	for intStimType=1:intStimTypes
		vecTrialsOfType = find(vecStimTypes==intStimType);
		matRespNSR(:,intStimType,1:numel(vecTrialsOfType)) = matResp(:,vecTrialsOfType);
	end
	
end
