function sTypes = getStimulusTypes(sIn,cellFields)
	%getStimulusTypes Retrieve unique stimulus combinations from data structure
	%   Syntax: sTypes = getStimulusTypes(sIn,[cellFields])
	%	Input can either be in ses, ses.structStim or sParTP format
	%
	%	Outputs three fields in structOut:
	%		matTypes		matrix containing all unique stimulus combinations
	%		vecNumTypes		vector containing number of unique stimuli per
	%						stimulus property
	%		cellNames		cell containing names of stimulus properties
	%
	%	Version history:
	%	1.0 - July 22 2013
	%	Created by Jorrit Montijn
	
	%% check which fields are present
	%field list
	if ~exist('cellFields','var')
		cellFields{1} = 'Orientation';
		cellFields{2} = 'Contrast';
		cellFields{3} = 'SpatialFrequency';
		cellFields{4} = 'Pitch';
		cellFields{5} = 'Scene';
		cellFields{6} = 'TemporalFrequency';
	end
	
	%get data location
	if isfield(sIn,'structStim') && (isfield(sIn.structStim,cellFields{1}) || isfield(sIn.structStim,'FrameOn'))
		sData = sIn.structStim;
	elseif isfield(sIn,cellFields{1})
		sData = sIn;
	end
	
	%pre-allocate variables
	cellNames = {};
	intUseFields = 0;
	matTypes = [];
	
	%loop through fields to check for presence
	for intField=1:numel(cellFields)
		strField = cellFields{intField};
		if isfield(sData,strField) && min(sData.(strField)) ~= max(sData.(strField))
			cellNames{end+1} = strField; %#ok<AGROW>
			intUseFields = intUseFields + 1;
		end
	end
	
	%% loop through fields to retrieve types
	vecNumTypes = zeros(1,intUseFields);
	for intField=1:intUseFields
		strField = cellNames{intField};
		
		%get stim values of all trials
		vecVals = sData.(strField);
		
		%retrieve uniques
		vecUniqueVals = unique(vecVals);
		
		%assign number to output vector
		intVals = length(vecUniqueVals);
		vecNumTypes(intField) = intVals;
		
		%assign uniques to indexing matrix
		if intField == 1
			matTypes = vecUniqueVals;
		else
			intMatLength = size(matTypes,2);
			matTypes = repmat(matTypes,[1 intVals]);
			vecAdd = sort(repmat(vecUniqueVals,[1 intMatLength]),'ascend');
			matTypes(intField,:) = vecAdd;
		end
	end
	
	%% assign output
	sTypes.matTypes = matTypes;
	sTypes.vecNumTypes = vecNumTypes;
	sTypes.cellNames = cellNames;
end

