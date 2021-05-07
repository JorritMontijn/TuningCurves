function [vecSpikeTimes,vecEventStarts] = getGeneratedTriPhasicR(m,Tr,T,tau,L_b,L_s)
	%getGeneratedTriPhasicR Generates triphasic neural data exponential ISI times
	%    [vecSpikeTimes,vecEventStarts] = getGeneratedTriPhasicR(m,Tr,T,tau,L_b,L_s)
	%
	%
	
	n_r = L_s*T;
	L_r = n_r/Tr;
	
	
	vecISIb = exprnd(1/L_b,[1 round(L_b*m*tau*5)]);
	vecBaseSpikeTimes = cumsum(vecISIb);
	vecISIs = exprnd(1/L_r,[1 round(L_r*m*tau*5)]);
	vecStimSpikeTimes = cumsum(vecISIs);
	cellStimSpikeTimesPerTrial = cell(1,m);
	cellBaseSpikeTimesPerTrial = cell(1,m);
	vecEventStarts = nan(m,1);
	for i_t=1:m
		tr_t0 = (i_t-1)*tau;
		tr_Tr = tr_t0+Tr;
		tr_T = tr_t0+T;
		tr_tau = tr_t0+tau;
		cellStimSpikeTimesPerTrial{i_t} = vecStimSpikeTimes(vecStimSpikeTimes>=tr_t0 & vecStimSpikeTimes<tr_Tr);
		cellBaseSpikeTimesPerTrial{i_t} = vecBaseSpikeTimes(vecBaseSpikeTimes>=tr_T & vecBaseSpikeTimes<tr_tau);
		vecEventStarts(i_t) = tr_t0;
	end
	vecSpikeTimes = sort(cat(1,cell2vec(cellStimSpikeTimesPerTrial),cell2vec(cellBaseSpikeTimesPerTrial)));
	
	
end

