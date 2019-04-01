function [dblDeltaPrime,dblSignificance,dblShuffledDeltaPrime,dblShuffledVariance] = testDeltaPrime(vecResp,vecAngles,intShuffleIters)
	%testDeltaPrime Calculates the significance of \(delta)'
	%	 [dblDeltaPrime,dblSignificance,dblShuffledVariance] = testDeltaPrime(vecResp,vecAngles,intShuffleIters)
	%
	%Inputs:
	% - vecResp; [1 x Trial] response vector
	% - vecAngles; [1 x Trial] stimulus orientation vector  [in radians]
	%
	%Note: you can also supply a matrix of responses [Neuron x Trial]
	%instead of a vector, and the function will return a vector of delta'
	%
	%Version History:
	%2019-03-13 Created delta-prime function [by Jorrit Montijn]
	
	%% check inputs
	if size(vecResp,2) ~= size(vecAngles,2)
		error([mfilename ':WrongInput'],'Response vector and angle vector are not the same length!');
	end
	if ~exist('intShuffleIters','var')
		intShuffleIters = 1000;
	end
	%%
	intC = intC+1;
	vecResp = vecResp(:,1:end-8);
	vecAngles = vecAngles(:,1:end-8);
	% calculate unshuffled
	dblDeltaPrime = getDeltaPrime(vecResp,vecAngles);
	
	% run shuffles
	vecDeltaPrimeShuffled = nan(size(vecResp,1),intShuffleIters);
	for intIter=1:intShuffleIters
		%% get shuffled data
		vecShuffledAngles = vecAngles(randperm(numel(vecAngles)));
		
		%% get delta prime
		dblDeltaPrimeShuffled = getDeltaPrime(vecResp,vecShuffledAngles);
	
		%% save to output
		vecDeltaPrimeShuffled(:,intIter) = dblDeltaPrimeShuffled;
	end
	
	% get variance
	dblShuffledDeltaPrime = mean(vecDeltaPrimeShuffled,2);
	dblShuffledVariance = std(vecDeltaPrimeShuffled,[],2);
	
	% get z-score
	dblVarDist = (dblDeltaPrime - dblShuffledDeltaPrime) ./ dblShuffledVariance;
	
	%get bias
	intN=size(vecResp,1);
	intReps=numel(vecAngles)/numel(unique(vecAngles));
	dblBC=(2*intReps-intN-3)/(2*intReps-2) - (2*intN)/intReps
	
	vecSize(intC) = intC;
	vecBias(intC) = (1-mean(dblShuffledDeltaPrime));
	vecBiasCorr(intC) = dblBC;
	vecCorrs(intC) = vecBias(intC) - vecBiasCorr(intC);
	
	
end

