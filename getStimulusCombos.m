function sCombos = getStimulusCombos(sIn,cellFields)
	%getStimulusCombos Retrieve unique stimulus combinations from data structure
	%   Syntax: sCombos = getStimulusCombos(sIn,cellFields)
	%	sIn must be structure with fields of same size, and values must be transformable to doubles
	%
	%Based on getStimulusTypes, works for >2 variables
	%
	%	Outputs three fields in sCombos, with n=trials, p=combos, q=variables
	%		vecUniqueVals	[1 x q] vector with # of unique values per variable
	%		matComboIdx		[p x q] matrix with index per variable for each combination
	%		matComboVal		[p x q] matrix with values per variable for each combination
	%		vecComboCounts	[p x 1] vector with # of trials for each combination
	%		matComboTrials	[p x n] matrix with included trials (boolean) for each combination
	%
	%	Version history:
	%	1.0 - November 9 2023
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
	
	%loop through fields to check for presence
	intUseFields = 0;
	cellNames = {};
	cellUniqueValsPerField = {};
	vecEntryNum = [];
	for intField=1:numel(cellFields)
		strField = cellFields{intField};
		if isfield(sData,strField) && min(sData.(strField)) ~= max(sData.(strField))
			intUseFields = intUseFields + 1;
			cellNames{intUseFields} = strField; %#ok<AGROW>
			cellUniqueValsPerField{intUseFields} = unique(sData.(strField));
			vecEntryNum(intUseFields) = numel(sData.(strField));
		end
	end
	if ~all(vecEntryNum(1)==vecEntryNum)
		error([mfilename ':InconsistentTrialNumbers'],'Not all fields have the same number of entries');
	end
	
	%% assign unique labels
	intTrialNum = vecEntryNum(1);
	vecUniqueVals = cellfun(@numel,cellUniqueValsPerField);
	intUniqueCombos = prod(vecUniqueVals);
	matComboIdx = nan(intUniqueCombos,intUseFields);
	matComboVal = nan(intUniqueCombos,intUseFields);
	vecComboCounts = zeros(intUniqueCombos,1);
	matComboTrials = false(intUniqueCombos,intTrialNum);
	
	%% loop through fields to retrieve types
	for intField=1:intUseFields
		if intField==1
			intFieldOffset=1;
		else
			intFieldOffset = prod(vecUniqueVals(1:(intField-1)));
		end
		matComboIdx(:,intField) = modx(ceil((1:intUniqueCombos)./intFieldOffset),vecUniqueVals(intField));
	end
	
	%% fill combo cells
	for intCombo=1:intUniqueCombos
		indTrials = true(intTrialNum,1);
		for intField=1:intUseFields
			matComboVal(intCombo,intField) = cellUniqueValsPerField{intField}(matComboIdx(intCombo,intField));
			strField = cellFields{intField};
			indTrials = indTrials & (flat(sData.(strField)) == matComboVal(intCombo,intField));
		end
		%add trial indices
		matComboTrials(intCombo,:) = indTrials;
		vecComboCounts(intCombo) = sum(indTrials);
	end
	
	%% assign output
	sCombos = struct;
	sCombos.vecUniqueVals = vecUniqueVals;
	sCombos.matComboIdx = matComboIdx;
	sCombos.matComboVal = matComboVal;
	sCombos.vecComboCounts = vecComboCounts;
	sCombos.matComboTrials = matComboTrials;
end

