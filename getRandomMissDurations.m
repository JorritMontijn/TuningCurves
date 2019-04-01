function sStimAggregate = getRandomMissDurations(sStimAggregate,intSwitch)
	%UNTITLED Summary of this function goes here
	%   Detailed explanation goes here
	
	if ~exist('intSwitch','var'),intSwitch = 1;end
	%calculate mean number of frames for detected stimuli
	intNumContrasts = length(unique(sStimAggregate.Contrast));
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
	
	%assign
	sStimAggregate.FrameOff = vecTrialDur + sStimAggregate.FrameOn;
end

