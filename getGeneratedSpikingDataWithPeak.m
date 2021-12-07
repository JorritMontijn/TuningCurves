function [vecSpikeTimes,dblPrefOri] = getGeneratedSpikingDataWithPeak(vecTrialAngles,matTrialT,dblBaseRate,dblPrefRate,dblJitter,dblKappa,boolDoublePeaked,dblPrefOri,intAddSpikes)
	%getGeneratedSpikingDataWithPeak Generates neural data using von Mises tuning curves
	%   [vecSpikeTimes,dblPrefOri] = getGeneratedSpikingDataWithPeak(vecTrialAngles,matTrialT,dblBaseRate,dblPrefRate,dblJitter,dblKappa,boolDoublePeaked,dblPrefOri,intAddSpikes)
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
	
	
	if ~exist('dblPrefOri','var') || isempty(dblPrefOri)
		dblPrefOri = rand(1)*2*pi;
	end
	if ~exist('boolDoublePeaked','var') || isempty(boolDoublePeaked)
		boolDoublePeaked = false;
	end
	if ~exist('intAddSpikes','var') || isempty(intAddSpikes)
		intAddSpikes = round(numel(matTrialT(:,1))/2);
	end
	
	%% get timings
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
	
	%generate peak response
	dblStartDelay = 0.1;
	vecTrialStarts = vecStarts(:)+dblStartDelay;
	vecChooseTrials = randperm(numel(vecStarts),intAddSpikes);
	%vecJitteredSpikes = (rand(size(vecChooseTrials))*0.002-0.001)*dblJitter;
	vecJitteredSpikes = 0.001*randn(size(vecChooseTrials))*dblJitter;
	vecJitteredSpikes = vecJitteredSpikes + mean(vecJitteredSpikes); %center spikes on 0
	vecPeakSpikes = vecJitteredSpikes(:) + vecTrialStarts(vecChooseTrials);
	
	%generate leading baseline
	vecITI = exprnd(1/dblBaseRate,[1 round(dblBaseRate*vecStarts(1)*10)]);
	vecSpT = cumsum(vecITI);
	vecLeadingSpikes = vecSpT(vecSpT < (vecStarts(1)+dblStartDelay));
	
	%generate responses by adding noise
	cellSpikeTimes = cell(intTrials,2);
	for intTrial=1:intTrials
		%stim
		vecITI = exprnd(1/vecMeanR(intTrial),[1 round(vecMeanR(intTrial)*vecStimDurs(intTrial)*10)]);
		vecSpT = cumsum(vecITI);
		cellSpikeTimes{intTrial,1} = vecSpT(vecSpT < vecStimDurs(intTrial))+vecStarts(intTrial)+dblStartDelay;
		
		%post-stim base base
		vecITI = exprnd(1/dblBaseRate,[1 round(dblBaseRate*vecBaseDurs(intTrial)*10)]);
		vecSpT = cumsum(vecITI);
		cellSpikeTimes{intTrial,2} = vecSpT(vecSpT < vecBaseDurs(intTrial))+vecStops(intTrial)+dblStartDelay;
	end
	%concatenate
	vecSpikeTimes = sort(cat(1,vecPeakSpikes(:),vecLeadingSpikes(:),cell2vec(cellSpikeTimes)),'ascend');
	
end

