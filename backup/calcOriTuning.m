function sTuning = calcOriTuning(ses,vecNeurons,structParams)
	%UNTITLED3 Summary of this function goes here
	%   Detailed explanation goes here
	%depends on:
	%ang2rad
	%circ_var
	%getStimulusTypes [depends on getUniqueVals.m]
	%getSelectionVectors
	%getNeuronResponse
	%imnorm
	
	%check input
	if ~exist('vecNeurons','var') || isempty(vecNeurons)
		vecNeurons = 1:numel(ses.neuron);
	end
	if nargin < 3
		structParams = struct;
	end
	
	%get orientation stimuli
	cellFields{1} = 'Orientation';
	sTypes = getStimulusTypes(ses,cellFields);
	cellSelect = getSelectionVectors(ses.structStim,sTypes);
	vecOrientations = sTypes.matTypes;
	
	%check if ori or dir
	if  (max(vecOrientations) - min(vecOrientations)) > 180
		boolOri = false;
		dblFactor = 1;
	else
		boolOri = true;
		dblFactor = 2;
	end
	
	% loop through stimulus types to get nr of reps
	intMaxReps = 0;
	for intS=1:length(cellSelect)
		intMaxReps = max(intMaxReps,sum(cellSelect{intS}));
	end
	
	%pre-allocate response matrix
	intOris = length(vecOrientations);
	intNeuronMax = max(vecNeurons);
	matStimResponse = nan(intOris,intMaxReps,intNeuronMax);
	
	%pre-allocate output
	vecOSI = nan(1,intNeuronMax);
	vecDSI = nan(1,intNeuronMax);
	vecPrefIndex = nan(1,intNeuronMax);
	vecPrefAngle = nan(1,intNeuronMax);
	
	% retrieve responses
	for intOriIndex=1:intOris
		vecSelect = cellSelect{intOriIndex};
		matThisStim = getNeuronResponse(ses, vecNeurons, vecSelect, structParams);
		intReps = size(matThisStim,2);
		if intReps < intMaxReps
			matThisStim = [matThisStim nan(length(vecNeurons),intMaxReps-intReps)];
		end
		matStimResponse(intOriIndex,:,:) = shiftdim(matThisStim,1);
	end
	
	%transform to ori (if dir)
	if boolOri
		matOriResponse = matStimResponse;
		vecOris = 2*ang2rad(vecOrientations);
	else
		intDirs = intOris;
		intOris = intDirs/2;
		vecFirstHalf = 1:intOris;
		vecSecondHalf = (intOris+1):intDirs;
		matOriResponse = (matStimResponse(vecFirstHalf,:,:) + matStimResponse(vecSecondHalf,:,:))/2;
		
		vecOris = 2*ang2rad(vecOrientations(vecFirstHalf));
	end
	
	%get tuning properties
	
	for intNeuron = vecNeurons
		%get data and normalize
		matOriResp = matOriResponse(:,:,intNeuron);
		vecOriResp = imnorm(mean(matOriResp,2));
		
		%get pref stim
		[dblPrefResp,intPrefIndex] = max(vecOriResp);
		dblPrefAngle = vecOrientations(intPrefIndex);
		vecPrefIndex(intNeuron) = intPrefIndex;
		vecPrefAngle(intNeuron) = dblPrefAngle;
		
		%get osi
		%dblCircVar = circ_var(vecOris', vecOriResp);
		%vecOSI(intNeuron) = 1 - dblCircVar;
		
		intOrthIndex = mod(intPrefIndex-1+round(intOris/2),intOris)+1;
		dblOrthResp = vecOriResp(intOrthIndex);
		vecOSI(intNeuron) = (dblPrefResp - dblOrthResp)/(dblPrefResp + dblOrthResp);
		
		%get dsi
		if ~boolOri
			%get data and normalize
			matTempResp = matStimResponse(:,:,intNeuron);
			vecTempResp = imnorm(nanmean(matTempResp,2));
			
			
			intOppIndex = mod(intPrefIndex-1+round(intDirs/2),intDirs)+1;
			dblOppResp = vecTempResp(intOppIndex);
			vecDSI(intNeuron) = 1 - abs(dblOppResp/dblPrefResp);
		end
	end
	
	%assign to output
	sTuning.vecOSI = vecOSI;
	sTuning.vecDSI = vecDSI;
	sTuning.vecPrefIndex = vecPrefIndex;
	sTuning.vecPrefAngle = vecPrefAngle;
end

