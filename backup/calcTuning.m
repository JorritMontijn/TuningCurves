function [structActivity,matAct_TypeByRep,vecStimTypeLookup,sOut] = calcTuning(ses,intNeuron,structParams)
	%calcTuning Calculates the tuning of a single neuron within a single
	%session.
	%syntax: [structActivity,matAct_TypeByRep] = calcTuning(ses,intNeuron,structParams)
	%	input:
	%	- ses: session structure
	%	- intNeuron: neuron to calculate tuning on
	%	- structParams: structure where you can set fields to define some
	%	parameters, including:
	%		-	structParams.doPlot
	%		-	structParams.startOffsetSecs
	%		-	structParams.vecStimIdentity
	%		-	structParams.vecPlotLookup
	%		-	structParams.vecStimTypeIndex
	%		-	structParams.vecStimTypeLookup
	%
	%	output:
	%	- structActivity: structure with dFoF matrix for each stimulus
	%	type and the corresponding frames, syntax:
	%		- structActivity(intStimTypeIndex).mat_dFoF
	%		- structActivity(intStimTypeIndex).vec_Frames(:,intRep) = [intStartFrame intStimOnFrame intStimOffFrame intStopFrame];
	%	- matAct_TypeByRep: 2D matrix of mean dFoF values per stim type and
	%	repetition number (intStimTypeIndex,intRep)
	%
	%Dependencies:
	% - none
	%
	%	Version history:
	%	1.0 - September 27 2012
	%	Created by Jorrit Montijn
	%
	%Note:
	%- ses structure should contain following fields:
	%	ses.structStim(intStim).Orientation
	%	ses.structStim(intStim).TrialNumber
	%	ses.structStim(intStim).FrameOn
	%	ses.structStim(intStim).FrameOff
	%	ses.neuron(intNeuron).dFoF
	
	%% inputs
	if ~exist('structParams','var')
		structParams = struct;
	end
	
	%defaults
	startOffsetSecs = -3; %before stim pres
	stopOffsetSecs = 2; %after stim pres
	vecStimIdentity = ses.structStim(:).Orientation;
	vecPlotLookup = [6 3 2 1 4 7 8 9];%
	vecIgnoreTrial = [];
	intRemoveRepetition = 0;
	
	vecStimTypeLookup = getOriListFromTrials(ses.structStim.Orientation);
	vecStimTypeIndex = 1:length(vecStimTypeLookup);
	
	%supplied
	if isfield(structParams,'doPlot'), doPlot = structParams.doPlot;else doPlot = false;end;
	if isfield(structParams,'startOffsetSecs'), startOffsetSecs = structParams.startOffsetSecs;end;
	if isfield(structParams,'vecStimIdentity'), vecStimIdentity = structParams.vecStimIdentity;end;
	if isfield(structParams,'vecPlotLookup'), vecPlotLookup = structParams.vecPlotLookup;end;
	if isfield(structParams,'vecStimTypeIndex'), vecStimTypeIndex = structParams.vecStimTypeIndex;end;
	if isfield(structParams,'vecStimTypeLookup'), vecStimTypeLookup = structParams.vecStimTypeLookup;end;
	if isfield(structParams,'vecIgnoreTrial')
		vecIgnoreTrial = structParams.vecIgnoreTrial;
		intRemoveRepetition = 2;
		if isfield(structParams,'boolRemoveRepetition'), intRemoveRepetition = structParams.boolRemoveRepetition;end;
	end
	
	%if isfield(structParams,'vecStimOris'), vecStimOris = structParams.vecStimOris;end;
	
	
	%make figure
	if doPlot,h=figure;end
	if ~isfield(structParams,'strType') || isempty(structParams.strType) || ~isfield(ses,structParams.strType)
		sCell = ses.neuron;
	else
		sCell = ses.(structParams.strType);
	end
	
	%calc offset in frames
	startOffset = ceil(startOffsetSecs*ses.samplingFreq);
	stopOffset = ceil(stopOffsetSecs*ses.samplingFreq);
	
	%pre-allocate
	intRepetitions = sum(ses.structStim.Orientation(1)==ses.structStim.Orientation);
	matAct_TypeByRep = nan(length(vecStimTypeIndex),intRepetitions);
	matStimAct_TypeByRep = nan(length(vecStimTypeIndex),intRepetitions);
	matBaseAct_TypeByRep = nan(length(vecStimTypeIndex),intRepetitions);
	
	%% splice by type
	for intStimTypeIndex=vecStimTypeIndex;
		%define vars
		intStimIdentity = vecStimTypeLookup(intStimTypeIndex);
		intOri = mod(intStimIdentity,360);
		vecStimLogical = vecStimIdentity == intStimIdentity;
		vecStim = find(vecStimLogical == 1);
		
		%figure
		if doPlot
			intPlotNumber = vecPlotLookup(intStimTypeIndex);
			subplot(3,3,intPlotNumber)
			if intStimTypeIndex == 1
				title(['Session=' ses.session '; Ori=' num2str(intOri) '; neuron=' num2str(intNeuron)]);
			else
				title(num2str(intOri));
			end
			hold on;
		end
		
		intReps = sum(vecStimLogical(:));
		for intRep=1:intReps
			intStim = vecStim(intRep);
			
			intStartFrame = ses.structStim.FrameOn(intStim) + startOffset;
			intStimOnFrame = ses.structStim.FrameOn(intStim);
			intStimOffFrame = ses.structStim.FrameOff(intStim);
			intStopFrame = ses.structStim.FrameOff(intStim) + stopOffset;
			
			vecActivity = sCell(intNeuron).dFoF(intStartFrame:intStopFrame);
			vecStimActivity = sCell(intNeuron).dFoF(intStimOnFrame:intStimOffFrame);
			
			structActivity(intStimTypeIndex).mat_dFoF(:,intRep) = vecActivity;
			structActivity(intStimTypeIndex).vec_Frames(:,intRep) = [intStartFrame intStimOnFrame intStimOffFrame intStopFrame];
			structActivity(intStimTypeIndex).mat_StimdFoF(:,intRep) = vecStimActivity;
			
			xVector = (startOffset/ses.samplingFreq):(1/ses.samplingFreq):((length(vecActivity)+startOffset-1)/ses.samplingFreq);
			
			%plot
			if doPlot
				plot(xVector,vecActivity,'k');
			end
			
			if ~ismember(intStim,vecIgnoreTrial)
				%output
				vecBaseAct = sCell(intNeuron).dFoF(intStartFrame:(intStimOnFrame-1));
				vecStimAct = sCell(intNeuron).dFoF(intStimOnFrame:intStimOffFrame);
				dblMeanAct = mean(vecStimAct) - mean(vecBaseAct);
				matAct_TypeByRep(intStimTypeIndex,intRep) = dblMeanAct;
				matStimAct_TypeByRep(intStimTypeIndex,intRep) = mean(vecStimAct);
				matBaseAct_TypeByRep(intStimTypeIndex,intRep) = mean(vecBaseAct);
			else
				%output
				vecBaseAct = nan(size(sCell(intNeuron).dFoF(intStartFrame:(intStimOnFrame-1))));
				vecStimAct = nan(size(sCell(intNeuron).dFoF(intStimOnFrame:intStimOffFrame)));
				dblMeanAct = nan(size(mean(vecStimAct) - mean(vecBaseAct)));
				matAct_TypeByRep(intStimTypeIndex,intRep) = nan(size(dblMeanAct));
				matStimAct_TypeByRep(intStimTypeIndex,intRep) = nan(size(mean(vecStimAct)));
				matBaseAct_TypeByRep(intStimTypeIndex,intRep) = nan(size(mean(vecBaseAct)));
			end
		end
		
		vecMean = nanmean(structActivity(intStimTypeIndex).mat_dFoF,2);
		
		
		if doPlot
			plot(xVector,vecMean,'b');
			xlim([min(xVector) max(xVector)]);
			xlim([-3 5]);
			ylim([-0.5 3]);
			hold off;
		end
	end
	if intRemoveRepetition == 1 %remove entire reptition if single trial is excluded
		[row,vecRem]=find(isnan(matAct_TypeByRep));
		vecCols = ~ismember(1:intReps,vecRem);
		matAct_TypeByRep = matAct_TypeByRep(:,vecCols);
		matStimAct_TypeByRep = matStimAct_TypeByRep(:,vecCols);
		matBaseAct_TypeByRep = matBaseAct_TypeByRep(:,vecCols);
	elseif intRemoveRepetition == 2 %remove random trials from stimulus classes (nr = max removed by user)
		matNan = isnan(matAct_TypeByRep);
		vecSum = sum(matNan,2);
		intRem = max(vecSum);
		vecRem = intRem - vecSum;
		if intRem == 1 %only max of 1; mathematically simpler
			for intStimType=find(vecRem>0)'
				intRemInd = randperm(intReps,1);
				matNan(intStimType,intRemInd) = true;
			end
		else %any number might have to be removed
			error
		end
		
		%remove nans
		matAct_TypeByRep(matNan) = nan;
		matAct_TypeByRep = sort(matAct_TypeByRep,2);
		matAct_TypeByRep = matAct_TypeByRep(:,1:(end-intRem));
		
		matStimAct_TypeByRep(matNan) = nan;
		matStimAct_TypeByRep = sort(matStimAct_TypeByRep,2);
		matStimAct_TypeByRep = matStimAct_TypeByRep(:,1:(end-intRem));
		
		matBaseAct_TypeByRep(matNan) = nan;
		matBaseAct_TypeByRep = sort(matBaseAct_TypeByRep,2);
		matBaseAct_TypeByRep = matBaseAct_TypeByRep(:,1:(end-intRem));
	end
	sOut.matStimAct_TypeByRep = matStimAct_TypeByRep;
	sOut.matBaseAct_TypeByRep = matBaseAct_TypeByRep;
	
	if doPlot
		subplot(3,3,5)
		theta = (vecStimTypeLookup*pi)/180;
		rho = mean(matAct_TypeByRep,2)';
		theta(end+1) = theta(1);
		rho(end+1) = rho(1);
		polar(theta,rho);
		set(gca,'XTick',[0:45:359])
	end
end