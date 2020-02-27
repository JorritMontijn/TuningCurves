function [vecR2,vecP] = doCrossValidatedTuningCurves(matData,vecStimOriDegrees,intTypeCV)
	%doCrossValidatedTuningCurves Compute tuning curve significance by cross-validated fitting
	%[dblR2,dblP] = ...
	%	doCrossValidatedTuningCurves(matData,vecTrialTypes,intTypeCV)
	%
	%Inputs:
	% - matData; [n x p]  Matrix of n observations/trials of p predictors/neurons
	% - vecTrialTypes; [n x 1] Trial indexing vector of c classes in n observations/trials
	% - intTypeCV; [int or vec] Integer switch 0-2 or trial repetition vector.
	%				Val=0, no CV; val=1, leave-one-out CV, val=2 (or
	%				vector), leave-repetition-out.
	%
	%Outputs:
	% - dblR2; [scalar] Fraction of correct classifications
	%
	%Version History:
	%2020-02-18 Created function [by Jorrit Montijn]
	
	%% check which kind of cross-validation
	if nargin < 3 || isempty(intTypeCV)
		intTypeCV = 2;
	end
	error('this is slow; try using means of stimuli to predict left out trial/rep')
	
	%% prepare
	intVerbose = 1;
	
	%get number of trials
	vecStimOriDegrees = vecStimOriDegrees(:);
	intTrials = numel(vecStimOriDegrees);
	
	%check if matData is [trial x neuron] or [neuron x trial]
	if size(matData,1) == intTrials && size(matData,2) == intTrials
		%number of neurons and trials is the same
		warning([mfilename ':SameNeuronsTrials'],'Number of neurons and trials is identical; please double check the proper orientation of [intNeurons x intTrials]');
	elseif size(matData,1) == intTrials
		%rotate
		matData = matData';
	elseif size(matData,2) == intTrials
		%size is correct
	else
		error([mfilename ':SameNeuronsTrials'],'Size of matData and vecTrialTypes do not match');
	end
	vecUniqueTrialTypes = unique(vecStimOriDegrees);
	intStimTypes = length(vecUniqueTrialTypes);
	vecTrialTypeIdx = label2idx(vecStimOriDegrees);
	
	%pre-allocate output
	ptrTic = tic;
	%vecStimOriDegrees =  vecTrialOris;
	
	%% check ori or dir
	intNeurons = size(matData,1);
	matData = double(matData);
	if range(vecStimOriDegrees) < (2*pi)
		warning([mfilename ':PossiblyRadians'],sprintf('Range of angles is %d degrees, are you sure you are not supplying radians?',range(vecStimOriDegrees)));
	elseif range(vecStimOriDegrees) > 180 %direction, full circle
		intParams = 5;
		vecKappa = [1 1];
		indKappa = logical([0 1 1 0 0]);
		funcFit = @vonMisesDoubleFitPX;
	else %orientation, half-circle
		intParams = 4;
		vecKappa = 1;
		indKappa = logical([0 1 0 0]);
		funcFit = @vonMisesSingleFitPX;
		vecStimOriDegrees = vecStimOriDegrees*2;
	end
	%get stimulus response by repetition; from [N x T] to [N x S x R]
	[matRespNSR,vecStimTypes,vecUniqueDegs] = getStimulusResponses(matData,vecStimOriDegrees);
	vecStimOriRads = ang2rad(vecStimOriDegrees);
	vecUniqueRads =  ang2rad(vecUniqueDegs);
	
	%% pre-allocate output
	sOptions = curvefitoptimoptions('curvefitfun','MaxFunEvals',1000,'MaxIter',1000,'Display','off');
	hTic = tic;
	
	%% cross-validate
	if numel(intTypeCV) == intTrials
		%third input is train/test set
		indTrainTrials = intTypeCV==min(intTypeCV);
		indTestTrials = find(intTypeCV==max(intTypeCV));
		matR2 = nan(intNeurons,1);
		
		for intNeuron=1:intNeurons
			%msg
			if toc(ptrTic) > 5
				ptrTic = tic;
				pause(eps);
				if intVerbose > 0,fprintf('Fitting; now at neuron %d/%d [%s]\n',intNeuron,intNeurons,getTime);end
			end
			
			%get data
			vecTrainData = matData(intNeuron,indTrainTrials);
			vecTrainOriRads = vecStimOriRads(indTrainTrials);
			vecTestData = matData(intNeuron,indTestTrials);
			vecTestOriRads = vecStimOriRads(indTestTrials);
			
			%prep initial vars
			[vecRespSR,vecStimTypes,vecUniqueDegs] = getStimulusResponses(vecTrainData,vecTrainOriRads);
			vecMeanRespPerOri = nanmean(vecRespSR(1,:,:),3);
			
			%build initial parameter vector
			dblPrefOri = mod(circ_mean(vecUniqueRads(:),vecMeanRespPerOri(:)),2*pi);
			dblBaseline = min(vecMeanRespPerOri);
			dblGain = range(vecMeanRespPerOri);
			vecP0 = [dblPrefOri vecKappa dblBaseline dblGain];
			
			%do fitting
			try
				vecFittedParams = lsqcurvefit(funcFit, vecP0, vecTrainOriRads, vecTrainData,[0 eps*vecKappa min(vecMeanRespPerOri) 0],[2*pi 1000*vecKappa mean(vecMeanRespPerOri) 1000],sOptions);
			catch
				vecFittedParams = curvefitfun(funcFit, vecP0, vecTrainOriRads, vecTrainData,[0 eps*vecKappa min(vecMeanRespPerOri) 0],[2*pi 1000*vecKappa mean(vecMeanRespPerOri) 1000],sOptions);
			end
			
			%calculate R^2
			vecFittedResp = feval(funcFit,vecFittedParams,vecTestOriRads);
			dblSSRes = sum((vecTestData - vecFittedResp).^2);
			dblSSTot = sum((vecTestData - mean(vecTrainData)).^2);
			matR2(intNeuron,1) = 1 - (dblSSRes / dblSSTot);
		end
		
	elseif intTypeCV == 0
		%no CV
		matR2 = nan(intNeurons,1);
		
		for intNeuron=1:intNeurons
			%msg
			if toc(ptrTic) > 5
				ptrTic = tic;
				pause(eps);
				if intVerbose > 0,fprintf('Fitting; now at neuron %d/%d [%s]\n',intNeuron,intNeurons,getTime);end
			end
			
			%prep initial vars
			vecMeanRespPerOri = nanmean(matRespNSR(intNeuron,:,:),3);
			%build initial parameter vector
			dblPrefOri = mod(circ_mean(vecUniqueRads(:),vecMeanRespPerOri(:)),2*pi);
			dblBaseline = min(vecMeanRespPerOri);
			dblGain = range(vecMeanRespPerOri);
			vecP0 = [dblPrefOri vecKappa dblBaseline dblGain];
			
			%do fitting
			try
				vecFittedParams = lsqcurvefit(funcFit, vecP0, vecStimOriRads, matData(intNeuron,:),[0 eps*vecKappa min(vecMeanRespPerOri) 0],[2*pi 1000*vecKappa mean(vecMeanRespPerOri) 1000],sOptions);
			catch
				vecFittedParams = curvefitfun(funcFit, vecP0, vecStimOriRads, matData(intNeuron,:),[0 eps*vecKappa min(vecMeanRespPerOri) 0],[2*pi 1000*vecKappa mean(vecMeanRespPerOri) 1000],sOptions);
			end
			
			%calculate R^2
			vecFittedResp = feval(funcFit,vecFittedParams,vecStimOriRads);
			dblSSRes = sum((matData(intNeuron,:) - vecFittedResp).^2);
			dblSSTot = sum((matData(intNeuron,:) - mean(matData(intNeuron,:))).^2);
			matR2(intNeuron,1) = 1 - (dblSSRes / dblSSTot);
		end
		
	elseif intTypeCV == 1
		%leave one out
		matR2 = nan(intNeurons,intTrials);
		
		for intNeuron=1:intNeurons
			%msg
			if toc(ptrTic) > 5
				ptrTic = tic;
				pause(eps);
				if intVerbose > 0,fprintf('Fitting; now at neuron %d/%d [%s]\n',intNeuron,intNeurons,getTime);end
			end
			
			%prep initial vars
			vecMeanRespPerOri = nanmean(matRespNSR(intNeuron,:,:),3);
			%build initial parameter vector
			dblPrefOri = mod(circ_mean(vecUniqueRads(:),vecMeanRespPerOri(:)),2*pi);
			dblBaseline = min(vecMeanRespPerOri);
			dblGain = range(vecMeanRespPerOri);
			vecP0 = [dblPrefOri vecKappa dblBaseline dblGain];
			
			
			%leave one out
			for intLeaveOut=1:intTrials
				%get info on to-be-left-out trial
				indTrain = ~isnan(vecStimOriDegrees);
				indTrain(intLeaveOut) = false;
				
				%split trials
				vecTrainData = matData(intNeuron,indTrain);
				vecTrainOriRads = vecStimOriRads(indTrain);
				vecTestData = matData(intNeuron,intLeaveOut);
				vecTestOriRads = vecStimOriRads(intLeaveOut);
				
				%do fitting
				try
					vecFittedParams = lsqcurvefit(funcFit, vecP0, vecTrainOriRads, vecTrainData,[0 eps*vecKappa min(vecMeanRespPerOri) 0],[2*pi 1000*vecKappa mean(vecMeanRespPerOri) 1000],sOptions);
				catch
					vecFittedParams = curvefitfun(funcFit, vecP0, vecTrainOriRads, vecTrainData,[0 eps*vecKappa min(vecMeanRespPerOri) 0],[2*pi 1000*vecKappa mean(vecMeanRespPerOri) 1000],sOptions);
				end
				
				%calculate R^2
				vecFittedResp = feval(funcFit,vecFittedParams,vecTestOriRads);
				dblSSRes = sum((vecTestData - vecFittedResp).^2);
				dblSSTot = sum((vecTestData - mean(vecTrainData)).^2);
				matR2(intNeuron,intLeaveOut) = 1 - (dblSSRes / dblSSTot);
			end
		end
		
	elseif intTypeCV == 2
		%cross validate by repetition
		[varDataOut,vecUnique,vecCounts,cellSelect,vecTrialRepetition] = label2idx(vecStimOriDegrees);
		intRepNum = max(vecTrialRepetition);
		matR2 = nan(intNeurons,intRepNum);
		
		for intNeuron=1:intNeurons
			%msg
			if toc(ptrTic) > 5
				ptrTic = tic;
				pause(eps);
				if intVerbose > 0,fprintf('Fitting; now at neuron %d/%d [%s]\n',intNeuron,intNeurons,getTime);end
			end
			
			%prep initial vars
			vecMeanRespPerOri = nanmean(matRespNSR(intNeuron,:,:),3);
			%build initial parameter vector
			dblPrefOri = mod(circ_mean(vecUniqueRads(:),vecMeanRespPerOri(:)),2*pi);
			dblBaseline = min(vecMeanRespPerOri);
			dblGain = range(vecMeanRespPerOri);
			vecP0 = [dblPrefOri vecKappa dblBaseline dblGain];
			
			%remove repetition
			for intRep=1:intRepNum
				
				%remove trials
				indThisRep = vecTrialRepetition==intRep;
				vecTrainData = matData(intNeuron,~indThisRep);
				vecTrainOriRads = vecStimOriRads(~indThisRep);
				vecTestData = matData(intNeuron,indThisRep);
				vecTestOriRads = vecStimOriRads(indThisRep);
				
				%do fitting
				try
					vecFittedParams = lsqcurvefit(funcFit, vecP0, vecTrainOriRads, vecTrainData,[0 eps*vecKappa min(vecMeanRespPerOri) 0],[2*pi 1000*vecKappa mean(vecMeanRespPerOri) 1000],sOptions);
				catch
					vecFittedParams = curvefitfun(funcFit, vecP0, vecTrainOriRads, vecTrainData,[0 eps*vecKappa min(vecMeanRespPerOri) 0],[2*pi 1000*vecKappa mean(vecMeanRespPerOri) 1000],sOptions);
				end
				
				%calculate R^2
				vecFittedResp = feval(funcFit,vecFittedParams,vecTestOriRads);
				dblSSRes = sum((vecTestData(:) - vecFittedResp(:)).^2);
				dblSSTot = sum((vecTestData(:) - mean(vecTrainData)).^2);
				matR2(intNeuron,intRep) = 1 - (dblSSRes / dblSSTot);
			end
		end
		
	else
		error([mfilename ':SyntaxError'],'CV type not recognized');
	end
	
	%output
	vecR2 = mean(matR2,2);
	[h,vecP] = ttest(matR2');
end

