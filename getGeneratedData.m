function [matResp,vecPrefOri] = getGeneratedData(intPopN,vecTrialAngles,vecKappa,dblDistDprime)
	%getGeneratedData Generates neural data using von Mises tuning curves
	%    matResp = getGeneratedData(intN,intRep,dblKappa,dblDistDprime)
	%
	%Inputs:
	% - intN; number of neurons to simulate
	% - intRep; number of repetitions to simulate
	% - vecTrialAngles; vector with angles at which to generate responses
	% - vecKappa; sharpness of tuning curve per neuron (von Mises Kappa)
	% - dblDistDprime; signal-to-noise in d' between pref and orth responses
	%
	%Outputs:
	% - matResp; [Neuron x Trial] response matrix with Gaussian noise
	% - vecTrialAngles; angles per trial
	% - vecPrefOri; preferred orientations used to generate responses
	%
	%Version History:
	%2019-03-14 Created getGeneratedData function [by Jorrit Montijn]
	
	%% generate preferred orientations
	vecPrefOri = rand(1,intPopN)*2*pi;
	intTrials = numel(vecTrialAngles);
	matResp = nan(intPopN,intTrials);
	if isscalar(vecKappa),vecKappa=repmat(vecKappa,[1 intPopN]);end
	for intN=1:intPopN
		%get mean tuning curve
		vecMeanR = circ_vmpdf(vecTrialAngles, vecPrefOri(intPopN), vecKappa(intPopN));% +...
		%	circ_vmpdf(vecTrialAngles, vecPrefOri(intPopN)+deg2rad(180), vecKappa(intPopN));
		%normalize to [10 dblDistDprime+10]
		vecMeanR = vecMeanR-min(vecMeanR);
		vecMeanR = dblDistDprime*(vecMeanR./max(vecMeanR));
		vecMeanR = vecMeanR+2;
		
		%generate responses by adding noise
		matResp(intN,:) = poissrnd(vecMeanR);
	end
end

