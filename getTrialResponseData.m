function [matTrialResponse,cellSelectContrasts,vecTrialDur] = getTrialResponseData(ses,sStimAggregate,intSwitch)
	%UNTITLED Summary of this function goes here
	%   Detailed explanation goes here
	
	%transform structure-based data to raw dFoF matrix
	if ~exist('intSwitch','var'),intSwitch = 1;end
	intNumContrasts = length(unique(ses.structStim.Contrast));
	intFrames = length(ses.neuron(1).dFoF);
	intNeurons = numel(ses.neuron);
	matActivity = zeros(intNeurons,intFrames);
	for intNeuron=1:intNeurons
		matActivity(intNeuron,:) = ses.neuron(intNeuron).dFoF;
	end
	
	%check for stim aggregate; if present, then calculate mean number of
	%frames for detected stimuli
	vecTrialDur = ses.structStim.FrameOff - ses.structStim.FrameOn;
	if exist('sStimAggregate','var') && isstruct(sStimAggregate)
		vecTrialDur = sStimAggregate.FrameOff - sStimAggregate.FrameOn;
		cellFieldsC_SA{1} = 'Contrast';
		sTypesC_SA = getStimulusTypes(sStimAggregate,cellFieldsC_SA);
		cellSelectContrasts = getSelectionVectors(sStimAggregate,sTypesC_SA);
		
		%get responded trials
		vecResponded = logical(sStimAggregate.vecTrialResponse);
		for intContrastIndex=1:intNumContrasts
			%get target trials
			vecContrastTrials = cellSelectContrasts{intContrastIndex};
			vecRespTrials = vecResponded & vecContrastTrials;
			vecNoRespTrials = ~vecResponded & vecContrastTrials;
			if sum(vecRespTrials) == 0 || sum(vecNoRespTrials) == 0
				continue;
			end
			
			%assign duration to no-response trials
			if intSwitch == 1
				vecRTs = sStimAggregate.FrameOff(vecRespTrials)-sStimAggregate.FrameOn(vecRespTrials);
				intRespTrials = sum(vecRespTrials);
				intNoRespTrials = sum(vecNoRespTrials);
				vecRT_Ext = repmat(vecRTs,[1 ceil(intNoRespTrials / intRespTrials)]);
				vecRT_M = vecRT_Ext(randperm(intNoRespTrials));
				vecTrialDur(vecNoRespTrials) = vecRT_M;
			else
				dblMeanResp = round(mean(sStimAggregate.FrameOff(vecRespTrials)-sStimAggregate.FrameOn(vecRespTrials)));
				vecTrialDur(vecNoRespTrials) = dblMeanResp;
			end
		end
	else
		sStimAggregate = ses.structStim;
		cellFieldsC_SA{1} = 'Contrast';
		sTypesC_SA = getStimulusTypes(sStimAggregate,cellFieldsC_SA);
		cellSelectContrasts = getSelectionVectors(sStimAggregate,sTypesC_SA);
	end
	
	%calculate mean per trial
	intTrials = length(vecTrialDur);
	matTrialResponse = zeros(intNeurons,intTrials);
	for intTrial = 1:intTrials
		intFrameOn = sStimAggregate.FrameOn(intTrial);
		intDur = vecTrialDur(intTrial);
		intFrameOff = intDur + intFrameOn;
		matTrialResponse(:,intTrial) = mean(matActivity(:,intFrameOn:intFrameOff),2);
	end
end

