function [sOut] = getTuningCurves(matResp,vecStimOriDegrees)
	%getTuningCurves Get tuning curves for neurons using von Mises fit.
	%Syntax:
	%[sOut] = getTuningCurves(matResp,vecStimOriDegrees)
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
	
	%% check ori or dir
	if range(vecStimOriDegrees) > 180 %direction, full circle
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
	[matRespNSR,vecStimTypes,vecUniqueDegs] = getStimulusResponses(matResp,vecStimOriDegrees);
	vecStimOriRads = ang2rad(vecStimOriDegrees);
	vecUniqueRads =  ang2rad(vecUniqueDegs);
	
	%% pre-allocate output
	intNeurons = size(matResp,1);
	matFittedParams = nan(intNeurons,intParams);
	matMeanResp = nan(intNeurons,numel(vecUniqueRads));
	matSDResp = nan(intNeurons,numel(vecUniqueRads));
	matFittedResp = nan(intNeurons,numel(vecUniqueRads));
	matVariance = nan(intNeurons,numel(vecKappa));
	matBandwidth = nan(intNeurons,numel(vecKappa));
	sOptions = curvefitoptimoptions('curvefitfun','MaxFunEvals',1000,'MaxIter',1000,'Display','off');
	hTic = tic;
	
	%% run neuron loop
	for intNeuron = 1:intNeurons
		%msg
		if toc(hTic) > 5
			hTic = tic;
			fprintf('  getTuningCurves: %d/%d [%s]\n',intNeuron,intNeurons,getTime);
		end
		
		%get responses
		vecResp = matResp(intNeuron,:);
		vecMeanRespPerOri = xmean(matRespNSR(intNeuron,:,:),3);
		vecSDRespPerOri = xstd(matRespNSR(intNeuron,:,:),3);
		
		%build initial parameter vector
		dblPrefOri = mod(circ_mean(vecUniqueRads(:),vecMeanRespPerOri(:)),2*pi);
		dblBaseline = min(vecMeanRespPerOri);
		dblGain = range(vecMeanRespPerOri);
		vecP0 = [dblPrefOri vecKappa dblBaseline dblGain];
		
		%do fitting & save data
		if exist('curvefitfun.m','file')
			matFittedParams(intNeuron,:) = curvefitfun(funcFit, vecP0, vecStimOriRads, vecResp,[0 eps*vecKappa min(vecMeanRespPerOri) 0],[2*pi 1000*vecKappa mean(vecMeanRespPerOri) 1000],sOptions);
		else
			matFittedParams(intNeuron,:) = lsqcurvefit(funcFit, vecP0, vecStimOriRads, vecResp,[0 eps*vecKappa min(vecMeanRespPerOri) 0],[2*pi 1000*vecKappa mean(vecMeanRespPerOri) 1000],sOptions);
		end
		matFittedResp(intNeuron,:) = feval(funcFit,matFittedParams(intNeuron,:),vecUniqueRads);
		matVariance(intNeuron,:) = 1 - (besseli(1,matFittedParams(intNeuron,indKappa)) ./ besseli(0,matFittedParams(intNeuron,indKappa))); %var=1-I1(k)/I0(k)
		matBandwidth(intNeuron,:) = 2*acos(1-((1./matFittedParams(intNeuron,indKappa))*log(2)));%FWHM=2*arccos(1- [(1/kappa) * ln(2)] )
		matMeanResp(intNeuron,:) = vecMeanRespPerOri;
		matSDResp(intNeuron,:) = vecSDRespPerOri;
		
		if boolPlot
			cla;
			errorbar(vecUniqueRads,vecMeanRespPerOri,vecSDRespPerOri./sqrt(size(matRespNSR,3)));
			hold on
			plot(vecUniqueRads,matFittedResp(intNeuron,:));
			hold off
			drawnow
			pause
		end
	end
	
	%% build output matrix
	sOut = struct;
	sOut.vecStimTypes = vecStimTypes;
	sOut.vecUniqueDegs = vecUniqueDegs;
	sOut.vecUniqueRads = vecUniqueRads;
	sOut.matMeanResp = matMeanResp;
	sOut.matSDResp = matSDResp;
	sOut.indKappa = indKappa;
	sOut.matFittedParams = matFittedParams;
	sOut.matFittedResp = matFittedResp;
	sOut.matVariance = matVariance;
	sOut.matBandwidth = matBandwidth;
end

