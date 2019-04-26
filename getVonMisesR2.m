function [vecCVFitR2,sOut] = getVonMisesR2(matResp,vecStimOriDegrees)
	%getTuningCurves Get tuning curves for neurons using von Mises fit.
	%Syntax:
	%[vecCVFitR2,sOut] = getVonMisesR2(matResp,vecStimOriDegrees)
	%
	%Inputs:
	% - matResp; [Neuron x Trial] response matrix
	% - vecStimOriDegrees; [1 x Trial] stimulus orientation vector  [in degrees]
	%
	%Output structure has fields:
	%sOut.vecStimTypes = indexed values corresponding to orientations
	%sOut.vecUniqueDegs = orientation in degrees (possibly x2 to complete
	%	the unit circle) in order of stimulus types as indexed by vecStimTypes
	%sOut.vecUniqueRads = same as vecUniqueDegs, but transformed to radians
	%sOut.matMeanResp = [Neuron x Orientation] matrix of mean responses per
	%	orientation, following indexing order of vecStimTypes
	%sOut.matSDResp = [Neuron x Orientation] matrix of standard deviation
	%	of responses per orientation, as above
	%sOut.indKappa = logical vector indicating which parameters are kappas
	%sOut.matFittedParams = [Neuron x parameters] matrix with von Mises
	%	parameters after fitting (4 or 5 parameters for ori/dir respectively)
	%sOut.matFittedResp = [Neuron x Orientation] matrix of fitted responses
	%sOut.matVariance = [Neuron x (1 or 2)] variance per neuron per
	%	response peak (1 for orientation, 2 for direction)
	%sOut.matBandwidth = [Neuron x (1 or 2)] bandwidth (Full-width-at-half-
	%	maximum) per neuron per response peak [in radians]
	%
	%	By Jorrit Montijn (Alex Pouget lab), 22-02-18 (dd-mm-yy; Universite de Geneve)
	
	%% header
	boolPlot=false;
	matResp = double(matResp);
	%vecStimOriDegrees =  vecTrialOris;
	
	%% check whether deg or rad
	if range(vecStimOriDegrees) <= 2*pi
		warning([mfilename ':RadLikely'],'Supplied angles are likely in radians, please supply degrees!');
	end
	%% check ori or dir
	if range(vecStimOriDegrees) > 180 %direction, full circle
		intParams = 5;
		vecKappaAndDir = [1 0.01];
		indKappaAndDir = logical([0 1 1 0 0]);
		funcFit = @vonMisesDoubleFitPX;
	else %orientation, half-circle
		intParams = 4;
		vecKappaAndDir = 1;
		indKappaAndDir = logical([0 1 0 0]);
		funcFit = @vonMisesSingleFitPX;
		vecStimOriDegrees = vecStimOriDegrees*2;
	end
	%get trial index
	[vecStimTypes,vecUniqueDegs,vecCounts,cellSelect,vecRepetition] = label2idx(vecStimOriDegrees);
	%equalize repetitions
	intTrials = numel(vecStimTypes);
	intReps = min(vecCounts);
	indRemTrials = vecRepetition>intReps;
	matResp(:,indRemTrials) = [];
	vecStimOriDegrees(indRemTrials) = [];
	%recalc trial idx
	[vecStimTypes,vecUniqueDegs,vecCounts,cellSelect,vecRepetition] = label2idx(vecStimOriDegrees);
	vecStimOriRads = ang2rad(vecStimOriDegrees);
	vecUniqueRads =  ang2rad(vecUniqueDegs);
	intNumOris = numel(vecUniqueDegs);
	intTrials = numel(vecStimTypes);
	intReps = min(vecCounts);
	
	%% pre-allocate output
	intNeurons = size(matResp,1);
	vecCVFitR2 = nan(intNeurons,1);
	matFittedParams = nan(intNeurons,intParams);
	matMeanResp = nan(intNeurons,intNumOris);
	matSDResp = nan(intNeurons,intNumOris);
	matFittedResp = nan(intNeurons,intNumOris);
	matVariance = nan(intNeurons,1);
	matBandwidth = nan(intNeurons,1);
	%strAlgorithm = 'levenberg-marquardt';
	strAlgorithm = 'trust-region-reflective';
	sOptions = curvefitoptimoptions('lsqcurvefit','MaxFunEvals',1000,'MaxIter',1000,'Algorithm',strAlgorithm,'Display','off');
	hTic = tic;
	
	%% run neuron loop
	for intNeuron = 1:intNeurons
		%msg
		if toc(hTic) > 5
			hTic = tic;
			fprintf('  getTuningCurves: %d/%d [%s]\n',intNeuron,intNeurons,getTime);
		end
		
		%get means/sds per stim type for initial parameters
		vecMeanRespPerOri = nan(1,intNumOris);
		vecSDRespPerOri = nan(1,intNumOris);
		for intStimType=1:intNumOris
			indSelect = cellSelect{intStimType};
			vecMeanRespPerOri(intStimType) = mean(matResp(intNeuron,indSelect));
			vecSDRespPerOri(intStimType) = std(matResp(intNeuron,indSelect));
		end
		
		%build initial parameter vector
		dblPrefOri = mod(circ_mean(vecUniqueRads(:),vecMeanRespPerOri(:)),2*pi);
		dblBaseline = min(vecMeanRespPerOri);
		dblGain = range(vecMeanRespPerOri);
		vecP0 = [dblPrefOri vecKappaAndDir dblBaseline dblGain];
		%get responses
		vecResp = matResp(intNeuron,:);
		%get initial fit
		%vecParams0 = lsqcurvefit(funcFit, vecP0, vecStimOriRads, vecResp,[0 eps*vecKappa min(vecMeanRespPerOri) 0],[2*pi 1000*vecKappa mean(vecMeanRespPerOri) 1000],sOptions);
		vecLB = [0 0.1*vecKappaAndDir min(vecMeanRespPerOri) 0];
		vecUB = [2*pi 100*vecKappaAndDir mean(vecMeanRespPerOri) 1000];
		vecParams0 = lsqcurvefit(funcFit, vecP0, vecStimOriRads, vecResp,vecLB,vecUB,sOptions);
		%assign overall fits to matrices
		matFittedParams(intNeuron,:) = vecParams0;
		matFittedResp(intNeuron,:) = feval(funcFit,matFittedParams(intNeuron,:),vecUniqueRads);
		matVariance(intNeuron,:) = 1 - (besseli(1,matFittedParams(intNeuron,find(indKappaAndDir,1))) ./ besseli(0,matFittedParams(intNeuron,find(indKappaAndDir,1)))); %var=1-I1(k)/I0(k)
		matBandwidth(intNeuron,:) = 2*acos(1-((1./matFittedParams(intNeuron,find(indKappaAndDir,1)))*log(2)));%FWHM=2*arccos(1- [(1/kappa) * ln(2)] )
		matMeanResp(intNeuron,:) = vecMeanRespPerOri;
		matSDResp(intNeuron,:) = vecSDRespPerOri;
		
		%leave one repetition out
		vecFitError = nan(1,intTrials);
		vecBaseError = nan(1,intTrials);
		%sOptions = curvefitoptimoptions('lsqcurvefit','MaxFunEvals',1000,'MaxIter',1000,'Display','off');
	
		for intLeaveRep=1:intReps
			indTestTrials = vecRepetition == intLeaveRep;
			%get training trials
			vecTrainResp = vecResp(~indTestTrials);
			vecTrainType = vecStimTypes(~indTestTrials);
			vecTrainRad = vecStimOriDegrees(~indTestTrials);
			%get test trials
			vecTestResp = vecResp(indTestTrials);
			vecTestType = vecStimTypes(indTestTrials);
			vecTestRad = vecStimOriDegrees(indTestTrials);
			
			
			%do fitting & save data
			%vecFittedParams = lsqcurvefit(funcFit, vecParams0, vecTrainRad, vecTrainResp,[0 0.1*vecKappa min(vecMeanRespPerOri) 0],[2*pi 100*vecKappa mean(vecMeanRespPerOri) 1000],sOptions);
			vecFittedParams = lsqcurvefit(funcFit, vecParams0, vecTrainRad, vecTrainResp,vecLB,vecUB,sOptions);
			vecFittedResp= feval(funcFit,vecFittedParams,vecTestRad);
			%get error
			vecFitError(indTestTrials) = vecTestResp - vecFittedResp;
			vecBaseError(indTestTrials) = vecTestResp - mean(vecTrainResp);
			
		end
		
		%calc R^2
		dblSS_tot = nansum(vecBaseError.^2);
		dblSS_res = nansum(vecFitError.^2);
		vecCVFitR2(intNeuron) = 1 - (dblSS_res / dblSS_tot);
	end
	
	%% build output matrix
	sOut = struct;
	sOut.vecStimTypes = vecStimTypes;
	sOut.vecUniqueDegs = vecUniqueDegs;
	sOut.vecUniqueRads = vecUniqueRads;
	sOut.matMeanResp = matMeanResp;
	sOut.matSDResp = matSDResp;
	sOut.indKappa = indKappaAndDir;
	sOut.matFittedParams = matFittedParams;
	sOut.matFittedResp = matFittedResp;
	sOut.matVariance = matVariance;
	sOut.matBandwidth = matBandwidth;
end

