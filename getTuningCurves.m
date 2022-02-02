function [sOut] = getTuningCurves(matResp,vecStimOriDegrees,boolPlot)
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
	if ~exist('boolPlot','var') || isempty(boolPlot),boolPlot=false;end
	matResp = double(matResp);
	%vecStimOriDegrees =  vecTrialOris;
	
	%% check ori or dir
	if range(vecStimOriDegrees) < (2*pi)
		warning([mfilename ':PossiblyRadians'],sprintf('Range of angles is %d degrees, are you sure you are not supplying radians?',range(vecStimOriDegrees)));
	elseif range(vecStimOriDegrees) > 180 %direction, full circle
		boolDouble = false;
		intParams = 5;
		vecKappa = [1 1];
		indKappa = logical([0 1 1 0 0]);
		funcFit = @vonMisesDoubleFitPX;
	else %orientation, half-circle
		boolDouble = true;
		intParams = 4;
		vecKappa = 1;
		indKappa = logical([0 1 0 0]);
		funcFit = @vonMisesSingleFitPX;
		vecStimOriDegrees = vecStimOriDegrees*2;
	end
	vecStimOriDegrees = vecStimOriDegrees(:)';
	%get stimulus response by repetition; from [N x T] to [N x S x R]
	[matRespNSR,vecStimTypes,vecUniqueDegs] = getStimulusResponses(matResp,vecStimOriDegrees);
	vecStimOriRads = ang2rad(vecStimOriDegrees);
	vecUniqueRads =  ang2rad(vecUniqueDegs);
	
	%% pre-allocate output
	intStimTypes = size(matRespNSR,2);
	intNeurons = size(matResp,1);
	matFittedParams = nan(intNeurons,intParams);
	matMeanResp = nan(intNeurons,numel(vecUniqueRads));
	matSDResp = nan(intNeurons,numel(vecUniqueRads));
	matFittedResp = nan(intNeurons,numel(vecUniqueRads));
	matVariance = nan(intNeurons,numel(vecKappa));
	matBandwidth = nan(intNeurons,numel(vecKappa));
	vecOriAnova = nan(intNeurons,1);
	intTests = (intStimTypes*(intStimTypes-1))/2;
	vecFitR2 = nan(intNeurons,1);
	vecFitT = nan(intNeurons,1);
	vecFitP = nan(intNeurons,1);
		
	try
		sOptions = optimoptions('curvefitfun','MaxFunEvals',1000,'MaxIter',1000,'Display','off');
	catch
		sOptions = curvefitoptimoptions('curvefitfun','MaxFunEvals',1000,'MaxIter',1000,'Display','off');
	end
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
		vecMeanRespPerOri = nanmean(matRespNSR(intNeuron,:,:),3);
		vecSDRespPerOri = nanstd(matRespNSR(intNeuron,:,:),[],3);
		
		%anova over orientations
		matR = squeeze(matRespNSR(intNeuron,:,:));
		vecR = matR(:);
		matG = repmat((1:intStimTypes)',[1 size(matR,2)]);
		vecG = matG(:);
		indRem = isnan(vecR);
		vecR(indRem) = [];
		vecG(indRem) = [];
		dblP_Anova = anova1(vecR,vecG,'off');
		vecOriAnova(intNeuron) = dblP_Anova;
		
		
		%build initial parameter vector
		dblPrefOri = mod(circ_mean(vecUniqueRads(:),vecMeanRespPerOri(:)),2*pi);
		dblBaseline = min(vecMeanRespPerOri);
		dblGain = range(vecMeanRespPerOri);
		vecP0 = [dblPrefOri vecKappa dblBaseline dblGain];
		
		%do fitting & save data
		%try
			matFittedParams(intNeuron,:) = lsqcurvefit(funcFit, vecP0, vecStimOriRads, vecResp,[0 eps*vecKappa min(vecMeanRespPerOri) 0],[2*pi 1000*vecKappa mean(vecMeanRespPerOri) 1000],sOptions);
		%catch
		%	matFittedParams(intNeuron,:) = curvefitfun(funcFit, vecP0, vecStimOriRads, vecResp,[0 eps*vecKappa min(vecMeanRespPerOri) 0],[2*pi 1000*vecKappa mean(vecMeanRespPerOri) 1000],sOptions);
		%end
		
		%get R^2
		vecFitR = feval(funcFit,matFittedParams(intNeuron,:),vecStimOriRads);
		%[dblR2,dblSS_tot,dblSS_res] = getR2(vecResp,vecFitR);
		[dblR2,dblSS_tot,dblSS_res,dblT,dblP] = getR2(vecResp,vecFitR,numel(vecP0));
		vecFitR2(intNeuron) = dblR2;
		vecFitT(intNeuron) = dblT;
		vecFitP(intNeuron) = dblP;
		
		%save output
		matFittedResp(intNeuron,:) = feval(funcFit,matFittedParams(intNeuron,:),vecUniqueRads);
		matVariance(intNeuron,:) = 1 - (besseli(1,matFittedParams(intNeuron,indKappa)) ./ besseli(0,matFittedParams(intNeuron,indKappa))); %var=1-I1(k)/I0(k)
		matBandwidth(intNeuron,:) = 2*acos(1-((1./matFittedParams(intNeuron,indKappa))*log(2)));%FWHM=2*arccos(1- [(1/kappa) * ln(2)] )
		matMeanResp(intNeuron,:) = vecMeanRespPerOri;
		matSDResp(intNeuron,:) = vecSDRespPerOri;
		
		
		if boolPlot
			cla;
			errorbar(rad2deg(vecUniqueRads),vecMeanRespPerOri,vecSDRespPerOri./sqrt(size(matRespNSR,3)));
			hold on
			plot(rad2deg(vecUniqueRads),matFittedResp(intNeuron,:));
			hold off
			title(sprintf('Pref deg=%.1f, R^2=%.3f,t=%.3f,p=%.2e',rad2deg(matFittedParams(intNeuron,1)),dblR2,dblT,dblP));
			drawnow
			%pause
		end
	end
	
	%% build output matrix
	if boolDouble
		%halve again
		vecUniqueDegs = vecUniqueDegs./2;
		vecUniqueRads = vecUniqueRads./2;
		matVariance = matVariance./2;
		matBandwidth = matBandwidth./2;
		matFittedParams(:,indKappa) = matFittedParams(:,indKappa)./2;
		matFittedParams(:,1) = matFittedParams(:,1)./2;
	end
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
	sOut.funcFit = funcFit;
	sOut.vecOriAnova = vecOriAnova;
	sOut.vecFitR2 = vecFitR2;
	sOut.vecFitT = vecFitT;
	sOut.vecFitP = vecFitP;
	
end

