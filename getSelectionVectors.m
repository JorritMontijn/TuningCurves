function cellSelect = getSelectionVectors(structStim,sTypes)
	%getSelectionVectors Retrieve indexing vectors
	%   Syntax: cellSelect = getSelectionVectors(structStim,sTypes)
	%	Uses output structure from getStimulusTypes() function as input
	%	(sTypes) as well as the structStim field from the ses structure.
	%	Outputs a cell array where each cell contains an indexing vector to
	%	select all trials/presentations beloning to each unique stimulus
	%	combination
	%
	%Dependencies:
	% - none
	%
	%	Version history:
	%	1.0 - July 22 2013
	%	Created by Jorrit Montijn
	%	1.1 - Feb 18 2020
	%	Trial selection from sTypes [by JM]
	
	%retrieve data from input structure
	matTypes = sTypes.matTypes;
	cellNames = sTypes.cellNames;

	%get trials
	intPresentations = length(structStim.(cellNames{1}));
	indBase = true(1,intPresentations);
		
	%loop through unique stimuli to create indexing vectors
	[intProperties,intStimuli] = size(matTypes);
	cellSelect = cell(1,intStimuli);
	for intCombo=1:intStimuli
		%select trials/presentations beloning to this stim combo
		indSelect = indBase;
		vecStimProps = matTypes(:,intCombo);
		
		%loop through stimulus modalities/properties
		for intProp=1:intProperties
			strName = cellNames{intProp};
			dblProp = vecStimProps(intProp);
			indProp = structStim.(strName) == dblProp;
			indSelect = indSelect & indProp;
		end
		
		%output data
		cellSelect{intCombo} = indSelect;
	end
end

