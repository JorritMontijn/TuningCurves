function matTuning = getTuningCurve( ses, vecNeurons )
    %Get tuning curves for neurons
    %
    %INPUT
    %- ses: session file from preprocessing
    %- vecNeurons: vector with neuron indexes (optional)
    %
    %OUTPUT
    %- matTuning: an NxRxO matrix of N neurons, with R repetitions and O
    %             orientations
    
    if nargin < 2
        vecNeurons = 1:numel(ses.neuron);
    end
    
    %Get trial selection vectors
    cellFields{1} = 'Orientation';
	sTypes = getStimulusTypes(ses,cellFields);
    sTypes.matTypes( sTypes.matTypes == 999 ) = [];
    sTypes.vecNumTypes = length( sTypes.matTypes );
	cellSelect = getSelectionVectors(ses.structStim, sTypes);
    
	% loop through stimulus types to get nr of reps
	intMaxReps = 0;
	for intS=1:length(cellSelect)
		intMaxReps = max(intMaxReps,sum(cellSelect{intS}));
	end
	
    %Loop through stim types and get response
    matTuning = nan(length(vecNeurons),intMaxReps,length(cellSelect));
    for s = 1:length(cellSelect)
        sParams = [];
%         sParams.intPreBaselineRemoval = 2;
		matThisStim = getNeuronResponse( ses, vecNeurons, cellSelect{s}, sParams );
		intReps = size(matThisStim,2);
		if intReps < intMaxReps
			matThisStim = [matThisStim nan(length(vecNeurons),intMaxReps-intReps)];
		end
        matTuning(:,:,s) = matThisStim;
    end
    matTuning = matTuning(vecNeurons,:,:);
    
end

