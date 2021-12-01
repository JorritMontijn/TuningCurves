function [vecSpikeTimes,dblPrefOri] = getGeneratedBurstingData(vecTrialRadians,matTrialT,sSpikingParams,sTuningParams)
	%getGeneratedBurstingData Generates neural data using von Mises tuning curves
	%   [vecSpikeTimes,dblPrefOri] = getGeneratedBurstingData(vecTrialAngles,matTrialT,sSpikingParams,sTuningParams)
	%
	%
	%Version History:
	%2021-03-09 Created getGeneratedBurstingData function [by Jorrit Montijn]
	%Parameters	ts (ms)	Number of bursts	Mean burst duration (ms)	Spikes in bursts/total spikes	Mean ISI in bursts (ms)
	%Typical bursts	34	17	104	622/625	2.91

	%% check inputs
	if ~exist('matTrialT','var') || isempty(matTrialT)
		if ~exist('vecTrialAngles','var') || isempty(vecTrialRadians)
			vecUniqueTrialAngles = 0:(360/24):359;
			vecTrialRadians = [];
			for intRep=1:10
				vecTrialRadians = cat(2,vecTrialRadians,vecUniqueTrialAngles(randperm(numel(vecUniqueTrialAngles))));
			end
			vecTrialRadians = deg2rad(vecTrialRadians);
		end
		intT = numel(vecTrialRadians);
		dblDur = 1;
		dblITI = 0.5;
		dblPreLead=1;
		vecStarts = dblPreLead:(dblDur+dblITI):((intT-0.5)*(dblDur+dblITI)+dblPreLead);
		vecStops = vecStarts + dblDur;
		matTrialT = cat(2,vecStarts',vecStops');
	end
	if ~exist('sSpikingParams','var') || isempty(sSpikingParams)
		sSpikingParams = [];
	end
	if ~exist('sTuningParams','var') || isempty(sTuningParams)
		sTuningParams = [];
	end
	
	%% spiking params
	dblBaseRate = getOr(sSpikingParams,'dblBaseRate',exprnd(1)); %mean baseline single spike rate (Hz) (exponential ISI)
	dblBurstEventRate = getOr(sSpikingParams,'dblBurstEventRate',exprnd(0.2)); %mean baseline rate of burst events (Hz) (exponential inter-event times)
	dblBurstDuration = getOr(sSpikingParams,'dblBurstDuration',normrnd(50,10)); %mean duration of burst events (ms) (Gamma, theta=0.5;k=2*dblBurstDuration)
	dblBurstISI = getOr(sSpikingParams,'dblBurstISI',(0.5+exprnd(2.9))/1000); %mean ISI during bursts (s) (Gamma, theta=0.5;k=2*dblBurstISI)
	
	%% tuning params
	boolDoublePeaked = getOr(sTuningParams,'boolDoublePeaked',false); %orientation or direction tuned
	dblPrefOri = getOr(sTuningParams,'dblPrefOri',rand(1)*2*pi); %preferred orientation (rads)
	dblKappa = getOr(sTuningParams,'dblKappa',rand(1)*5+5); %von Mises concentration parameter
	dblPrefRate = getOr(sTuningParams,'dblPrefRate',dblBaseRate); %mean single-spike rate during stimulus (exponential ISI)
	dblPrefBurstEventRate = getOr(sTuningParams,'dblPrefBurstEventRate',dblBurstEventRate+exprnd(2)); %mean evoked rate of burst events (Hz) (exponential inter-event times)
	
	%% inputs
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
	if ~exist('dblPrefOri','var') || isempty(dblPrefOri)
		dblPrefOri = rand(1)*2*pi;
	end
	intTrials = numel(vecTrialRadians);
	
	%get mean tuning curve
	vecMeanR = circ_vmpdf(vecTrialRadians, dblPrefOri, dblKappa);
	if boolDoublePeaked
		vecMeanR = vecMeanR + circ_vmpdf(vecTrialRadians, dblPrefOri+deg2rad(180), dblKappa);
	end
	
	%bursts
	vecBurstR = vecMeanR;
	if range(vecBurstR) == 0
		vecBurstR(:) = dblPrefBurstEventRate;
	else
		vecBurstR = vecBurstR ./ max(vecBurstR);
		vecBurstR = vecBurstR*(dblPrefBurstEventRate-dblBurstEventRate) + dblBurstEventRate;
	end
	
	%single spikes
	if range(vecMeanR) == 0
		vecMeanR(:) = dblPrefRate;
	else
		vecMeanR = vecMeanR-min(vecMeanR);
		vecMeanR = vecMeanR ./ max(vecMeanR);
		vecMeanR = vecMeanR*(dblPrefRate-dblBaseRate) + dblBaseRate;
	end
	
	%% build single-spike behavior
	%generate leading baseline
	vecITI = exprnd(1/dblBaseRate,[1 round(dblBaseRate*vecStarts(1)*10)]);
	vecSpT = cumsum(vecITI);
	vecLeadingSpikes = vecSpT(vecSpT < (vecStarts(1)+0.05+0.02*rand(1)));
	
	%generate responses by adding noise
	cellSpikeTimes = cell(intTrials,2);
	for intTrial=1:intTrials
		%% single-spike event
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
	
	%% add burst event times
	dblEndT = vecStops(end)+0.05+0.02*rand(1);
	vecBaseBurstITIs = exprnd(1/dblBurstEventRate,[1 round(dblBurstEventRate*dblEndT*10)]);
	vecBaseBurstT = cumsum(vecBaseBurstITIs);
	vecBaseBurstT(vecBaseBurstT>dblEndT)=[];
	cellBurstTimes = cell(intTrials,1);
	
	for intTrial=1:intTrials
		%remove baseline
		%vecBaseBurstT(vecBaseBurstT>vecStarts(intTrial) & vecBaseBurstT<vecStops(intTrial)) = [];
		
		%stim
		vecIBI = exprnd(1/vecBurstR(intTrial),[1 round(vecBurstR(intTrial)*vecStimDurs(intTrial)*10)]);
		vecBurstT = cumsum(vecIBI);
		vecBurstT = vecBurstT(vecBurstT < vecStimDurs(intTrial))+vecStarts(intTrial)+0.05+0.02*rand(1);
		cellBurstTimes{intTrial} = vecBurstT;
	end
	vecAllBursts = cat(1,vecBaseBurstT',cell2vec(cellBurstTimes));
	theta = 0.5;
	k_dur = 2*dblBurstDuration;
	vecBurstDur = gamrnd(k_dur,theta,size(vecAllBursts))/1000;
	
	%% remove overlapping bursts
	boolDone = false;
	indKeepBursts = true(size(vecAllBursts));
	intLastBurst = 1;
	while ~boolDone
		dblLastBurstStart = vecAllBursts(intLastBurst);
		dblLastBurstDur = vecBurstDur(intLastBurst);
		intNextBurst = find(vecAllBursts>dblLastBurstStart+dblLastBurstDur,1);
		if isempty(intNextBurst),break;end
		indKeepBursts((intLastBurst+1):(intNextBurst-1))=false;
		intLastBurst = intNextBurst;
	end
	vecAllBursts = vecAllBursts(indKeepBursts);
	vecBurstDur = vecBurstDur(indKeepBursts);
	
	%% build spike trains per burst
	intBursts = numel(vecAllBursts);
	cellBurstSpikeT = cell(intBursts,1);
	for intBurst=1:intBursts
		%dblBurstDuration = getOr(sSpikingParams,'dblBurstDuration',exprnd(0.1)); %mean duration of burst events (s) (Gamma, theta=0.5;k=2*dblBurstDuration)
		%dblBurstISI = getOr(sSpikingParams,'dblBurstISI',exprnd(2.9)/1000); %mean ISI during bursts (s) (Gamma, theta=0.5;k=2*dblBurstISI)
		dblBurstStart = vecAllBursts(intBurst);
		
		dblBurstDur = vecBurstDur(intBurst);
		
		theta = 0.5;
		k_isi = 2*dblBurstISI*1000;
		vecSpikeIntervals = gamrnd(k_isi,theta,[1 round((dblBurstDur/(dblBurstISI))*10)]);
		vecSpikeT = cumsum([0 vecSpikeIntervals])/1000;
		cellBurstSpikeT{intBurst} = vecSpikeT(vecSpikeT < dblBurstDur)+dblBurstStart;
	end
	
	%% concatenate everything
	vecSpikeTimes = sort(cat(1,vecSpikeTimes,cell2vec(cellBurstSpikeT)),'ascend');
end

