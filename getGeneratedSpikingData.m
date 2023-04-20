function [vecSpikeTimes,dblPrefOri] = getGeneratedSpikingData(vecTrialAngles,matTrialT,dblBaseRate,dblPrefRate,dblKappa,boolDoublePeaked,dblPrefOri)
	%getGeneratedData Generates neural data using von Mises tuning curves
	%    matResp = getGeneratedSpikingData(intN,intRep,dblKappa,dblDistDprime)
	%
	%
	%Version History:
	%2019-09-23 Created getGeneratedSpikingData function [by Jorrit Montijn]
	
	%% inputs
	%vecTrialAngles
	%vecTrialStarts
	%dblBaseRate
	%dblPrefRate
	%dblKappa
	%boolDoublePeaked
	
	
	if ~exist('boolDoublePeaked','var') || isempty(boolDoublePeaked)
		boolDoublePeaked = false;
	end
	vecStarts = matTrialT(:,1);
	if size(matTrialT,2) == 2
		vecStops = matTrialT(:,2);
	else
		vecStops = vecStarts + median(diff(vecStarts));
	end
	vecStimDurs = vecStops-vecStarts;
	vecBaseDurs = vecStarts(2:end) -  vecStops(1:(end-1));
	vecBaseDurs(end+1) = median(vecBaseDurs);
	%% generate preferred orientations
	if range(vecTrialAngles) > (2*pi)
		error([mfilename ':AnglesNotInRadians'],'Range of angles is >2pi');
	end
	if ~exist('dblPrefOri','var') || isempty(dblPrefOri)
		dblPrefOri = rand(1)*2*pi;
	end
	intTrials = numel(vecTrialAngles);
	
	%get mean tuning curve
	vecMeanR = circ_vmpdf(vecTrialAngles, dblPrefOri, dblKappa);
	if boolDoublePeaked
		vecMeanR = vecMeanR + circ_vmpdf(vecTrialAngles, dblPrefOri+deg2rad(180), dblKappa);
	end
	
	%normalize
	if range(vecMeanR) == 0
		vecMeanR(:) = dblPrefRate;
	else
		vecMeanR = vecMeanR-min(vecMeanR);
		vecMeanR = vecMeanR ./ max(vecMeanR);
		vecMeanR = vecMeanR*(dblPrefRate-dblBaseRate) + dblBaseRate;
	end
	
	%generate leading baseline
	vecITI = exprnd(1/dblBaseRate,[1 round(dblBaseRate*vecStarts(1)*10)]);
	vecSpT = cumsum(vecITI);
	vecLeadingSpikes = vecSpT(vecSpT < (vecStarts(1)+0.05+0.02*rand(1)));
	
	%generate responses by adding noise
	cellSpikeTimes = cell(intTrials,2);
	for intTrial=1:intTrials
		%stim
		vecITI = exprnd(1/vecMeanR(intTrial),[1 round(vecMeanR(intTrial)*vecStimDurs(intTrial)*10)]);
		vecSpT = cumsum(vecITI);
		cellSpikeTimes{intTrial,1} = vecSpT(vecSpT < vecStimDurs(intTrial))+vecStarts(intTrial)+0.05+0.02*rand(1);
		
		%post-stim base base
		vecITI = exprnd(1/dblBaseRate,[1 round(dblBaseRate*vecBaseDurs(intTrial)*10)]);
		vecSpT = cumsum(vecITI);
		cellSpikeTimes{intTrial,2} = vecSpT(vecSpT < vecBaseDurs(intTrial))+vecStops(intTrial)+0.05+0.02*rand(1);
	end
	%concatenate
	vecSpikeTimes = cat(1,vecLeadingSpikes',cell2vec(cellSpikeTimes));
end

