function vecSpikeTimes = getGeneratedMultiPhasicR(dblBaseRate,vecRates,vecDurs,vecRepStarts,dblEndT)
	%getGeneratedMultiPhasicR Generates multiphasic neural data with exponential ISI times
	%    vecSpikeTimes = getGeneratedMultiPhasicR(dblBaseRate,vecRates,vecDurs,vecRepStarts,dblEndT)
	%
	%dblBaseRate: baseline rate, additive to vecRates
    %vecRates: increase in rate during different epochs
    %vecDurs: durations of said epochs
    %vecRepStarts: vector with start times of repeats
    %
    %Example:
    %dblBaseRate = 5
    %vecRates = [1 3 0 5]
    %vecDurs = [2 1 1.5 1]
    %vecRepStarts = [5 11 21]
    %dblEndT = 30
    %these parameters will generate the response of neuron that fires at 5 Hz from t=0-5, 6 Hz from
    %t=5-7, 8 Hz from t=7-8, 5 Hz from t=8-9.5, 10 Hz from t=9.5-10.5, 5 Hz from t=10.5-11, and then
    %repeats the cycle from t=11-16.5, baseline from t=16.5-21, and a third repetition from
    %t=21-26.5, and an end period from t=26.5-30.
    %
    %If periods are overlapping (e.g., when setting vecRepStarts to [1 2 3]), then rates are
    %additive in the sense that the expectation of the ISI during that epoch is 1/additive_rate
    

    %% generate base rate
    vecBaseSpikes = cumsum(exprnd(1/dblBaseRate,[1 round(dblBaseRate*dblEndT*5)]));
	vecBaseSpikes(vecBaseSpikes>dblEndT)=[];

    cellTrialEpochSpikes = cell(numel(vecRepStarts),numel(vecDurs));
	for intTrial=1:numel(vecRepStarts)
        dblStartEpoch = vecRepStarts(intTrial);
        for intEpoch=1:numel(vecDurs)
            dblDur = vecDurs(intEpoch);
            dblEpochRate = vecRates(intEpoch);
            vecEpochSpikes = cumsum(exprnd(1/dblEpochRate,[1 round(dblEpochRate*dblDur*5)]));
            vecEpochSpikes(vecEpochSpikes > dblDur) = [];
	        cellTrialEpochSpikes{intTrial,intEpoch} = dblStartEpoch+vecEpochSpikes;
            dblStartEpoch = dblStartEpoch + dblDur;
        end
		
	end
	vecSpikeTimes = sort(cat(1,vecBaseSpikes(:),cell2vec(cellTrialEpochSpikes)));
end

