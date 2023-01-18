function structOut = calcStimCorrsNatMov(ses,varargin)
	%calcStimCorrs Calculates signal+noise correlations from stimulus
	%presentations
	%syntax: structOut = calcStimCorrs(ses)
	%	input:
	%	- ses: structure containing session information (prepro output)
	%
	%	output:
	%	- structOut: a structure containing the following fields:
	%		- matSignalCorrs: matrix containing signal correlations
	%		- matNoiseCorrs: matrix containing noise correlations
	%		- matStimResponse: matrix containing mean response of neuron by stimulus presentation
	%		- matSignalResponse: matrix containing mean response of neuron by stimulus type
	%		- sTypes: structure containing data on stimulus types of the
	%			different separated classes
	%
	%Dependencies:
	% - getStimulusTypes.m
	% - getSelectionVectors.m
	% - getNeuronResponse.m
	%
	%	Version history:
	%	1.0 - May 31 2013
	%	Created by Jorrit Montijn
	%	2.0 - July 25 2013
	%	Modified to work with arbitrary stimulus categories based on
	%	getStimulusTypes/getSelectionVectors/getNeuronResponse triple combo
	
	%% get input
	if nargin > 1
		sParams = varargin{1};
	end
	if ~exist('sParams','var'), sParams = struct;end
	if isfield(sParams,'cellFields'), cellFields = sParams.cellFields;else cellFields{1} = 'Orientation';end
	if isfield(sParams,'structParams'), structParams = sParams.structParams;else structParams = struct;end
	
	%% get stim response matrix
	%split data into bins
	intFramesPerBin = 12;
	vecOn = ses.structStim.FrameOn;
	vecOff = ses.structStim.FrameOff;
	vecDur = vecOff-vecOn;
	intMinDur = min(vecDur);
	vecBinOff = intFramesPerBin:intFramesPerBin:intMinDur;
	vecBinOn = vecBinOff-intFramesPerBin+1;
	
	%% pre-allocate response matrix
	intTypes = numel(vecBinOn);
	intReps = numel(vecOn);
	intNeurons = numel(ses.neuron);
	matStimResponse = nan(intTypes,intReps,intNeurons);
	
	%% retrieve responses
	for intNeuron=1:intNeurons
		vecdFoF=ses.neuron(intNeuron).dFoF;
		for intType=1:intTypes
			vecBinsOn = vecOn - 1.5 + vecBinOn(intType);
			vecBinsOff = vecOn - 0.5 + vecBinOff(intType);
			[vecCounts,vecMeans] = makeBins(1:numel(vecdFoF),vecdFoF,sort([vecBinsOn vecBinsOff]));
			matStimResponse(intType,:,intNeuron) = vecMeans(1:2:end);
		end
	end
	matStimResponse(:,1,:) = [];
	
	%% calculate signal correlations
	matSignalResponse = squeeze(mean(matStimResponse,2));%nan(intStims,intNeurons);
	
	%get mean per stim type
	matSignalCorrs = corr(matSignalResponse);
	
	%set autocorrelation at 0
	matSignalCorrs(diag(true(1,intNeurons))) = 0;
	
	
	%% do same for noise correlations
	matAllNoise = nan(intNeurons,intNeurons,intTypes);
	for intType=1:intTypes
		%calculate noise correlations per stimulus type
		matNoiseResponse = squeeze(matStimResponse(intType,:,:));
		matAllNoise(:,:,intType) = corr(matNoiseResponse);
	end
	%then average
	matNoiseCorrs = mean(matAllNoise,3);
	
	%set autocorrelation to 0
	matNoiseCorrs(diag(true(1,intNeurons))) = 0;
	
	%% put in output
	structOut.matAllNoise = matAllNoise;
	structOut.matSignalCorrs = matSignalCorrs;
	structOut.matNoiseCorrs = matNoiseCorrs;
	structOut.matStimResponse = matStimResponse;
	structOut.matSignalResponse = matSignalResponse;
end
