function structOut = calcStimCorrsAsym(ses,varargin)
	%calcStimCorrsAsym Calculates signal+noise correlations from stimulus
	%presentations, also works when number of repetitions per type are
	%inconsistent; it weights the noise correlation values obtained by the
	%number of observations they are are based on
	%syntax: structOut = calcStimCorrsAsym(ses)
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
	%	1.0 - August 29 2014
	%	Created by Jorrit Montijn, modified version of calcStimCorrs()
	
	%% get input
	if nargin > 1
		sParams = varargin{1};
	end
	if ~exist('sParams','var'), sParams = struct;end
	if isfield(sParams,'cellFields'), cellFields = sParams.cellFields;else cellFields = {'Orientation'};end
	if isfield(sParams,'structParams'), structParams = sParams.structParams;else structParams = struct;end
	
	%% get indexing vectors for unique stimulus combinations
	sTypes = getStimulusTypes(ses,cellFields);
	cellSelect = getSelectionVectors(ses.structStim,sTypes);
	intTypes = length(cellSelect);
	vecKeep = true(1,intTypes);
	for intType = 1:intTypes
		if sum(cellSelect{intType}) == 0
			vecKeep(intType) = false;
		end
	end
	cellSelect = cellSelect(vecKeep);
	sTypes.matTypes = sTypes.matTypes(:,vecKeep);
	
	%% pre-allocate response matrix
	intTypes = sum(vecKeep);
	intReps = sum(cellSelect{1});
	intNeurons = numel(ses.neuron);
	matStimResponse = nan(intTypes,intReps,intNeurons);
	
	%% retrieve responses
	for intType=1:intTypes
		vecSelect = cellSelect{intType};
		if sum(vecSelect) > 0
			matStimResponse(intType,1:sum(vecSelect),:) = shiftdim(getNeuronResponse(ses,1:intNeurons,vecSelect,structParams),1);
		end
	end
	matStimResponse(matStimResponse==0) = nan;
	
	%% calculate signal correlations
	matSignalResponse = squeeze(nanmean(matStimResponse,2));%nan(intStims,intNeurons);
	
	%get mean per stim type
	matSignalCorrs = corr(matSignalResponse);
	
	%set autocorrelation at 0
	matSignalCorrs(diag(true(1,intNeurons))) = 0;
	
	
	%% do same for noise correlations
	matAllNoise = nan(intNeurons,intNeurons,intTypes);
	vecWeighting = nan(1,intTypes);
	for intType=1:intTypes
		%calculate noise correlations per stimulus type
		matNoiseResponse = squeeze(matStimResponse(intType,:,:));
		indNans = isnan(matNoiseResponse(:,1));
		if sum(~indNans) < 2
			vecWeighting(intType) = 0;
			continue;
		end
		vecWeighting(intType) = sum(~indNans)-1;
		matNoiseResponse = matNoiseResponse(~indNans,:);
		matAllNoise(:,:,intType) = corr(matNoiseResponse)*vecWeighting(intType);
	end
	%then average
	matNoiseCorrs = nansum(matAllNoise,3)/sum(vecWeighting);
	
	%set autocorrelation to 0
	matNoiseCorrs(diag(true(1,intNeurons))) = 0;
	
	%% put in output
	structOut.matSignalCorrs = matSignalCorrs;
	structOut.matNoiseCorrs = matNoiseCorrs;
	structOut.matStimResponse = matStimResponse;
	structOut.matSignalResponse = matSignalResponse;
	structOut.sTypes = sTypes;
	structOut.cellSelect = cellSelect;
end
