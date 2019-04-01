function sOut = getTuningProperties(ses)
	%UNTITLED4 Summary of this function goes here
	%   Detailed explanation goes here
	
	%def variables
	intNeurons = numel(ses.neuron);
	matSelect = tril(true(intNeurons,intNeurons),-1);
	
	%distance
	sParams.vecEpochStart = 1;
	sParams.intEpochDuration = 2;
	
	%get signal+noise corrs
	structStimCorrs = calcStimCorrs(ses);
	
	%OSI, DSI, pref dir
	sTuning = calcOriTuning(ses);
	
	%get matrices
	matSignalCorrs = structStimCorrs.matSignalCorrs;
	matNoiseCorrs = structStimCorrs.matNoiseCorrs;
	sParams.boolAngDiff = false;
	sParams.boolMean = false;
	matDiffOSI = vec2diffmat(sTuning.vecOSI,sParams);
	matDiffDSI = vec2diffmat(sTuning.vecDSI,sParams);
	sParams.boolAngDiff = true;
	matDiffPrefAngle = rad2ang(vec2diffmat(ang2rad(sTuning.vecPrefAngle),sParams));
	sParams.boolMean = true;
	sParams.boolAngDiff = false;
	matMeanOSI = vec2diffmat(sTuning.vecOSI,sParams);
	matMeanDSI = vec2diffmat(sTuning.vecDSI,sParams);
	
	%transform to vectors
	sOut.vecDistributionSC = matSignalCorrs(matSelect);
	sOut.vecDistributionNC = matNoiseCorrs(matSelect);
	sOut.vecDistribution_dOSI = abs(matDiffOSI(matSelect));
	sOut.vecDistribution_dDSI = abs(matDiffDSI(matSelect));
	sOut.vecDistribution_dPA = abs(matDiffPrefAngle(matSelect));
	sOut.vecDistribution_mOSI = matMeanOSI(matSelect);
	sOut.vecDistribution_mDSI = matMeanDSI(matSelect);
end

