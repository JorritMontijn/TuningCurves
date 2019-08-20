function [matResp,vecPrefOri] = getGeneratedData(intPopN,vecTrialAngles,vecKappa,dblDistDprime,dblGaussNoise,boolDoublePeaked)
	%getGeneratedData Generates neural data using von Mises tuning curves
	%    matResp = getGeneratedData(intN,intRep,dblKappa,dblDistDprime)
	%
	%Inputs:
	% - intN; number of neurons to simulate
	% - intRep; number of repetitions to simulate
	% - vecTrialAngles; vector with angles at which to generate responses
	% - vecKappa; sharpness of tuning curve per neuron (von Mises Kappa)
	% - dblDistDprime; signal-to-noise in d' between pref and orth responses
	% - dblGaussNoise; additive mean-zero Gaussian noise (sd)
	% - boolDoublePeaked; orientation or direction tuned?
	%
	%Outputs:
	% - matResp; [Neuron x Trial] response matrix with Gaussian noise
	% - vecTrialAngles; angles per trial
	% - vecPrefOri; preferred orientations used to generate responses
	%
	%Version History:
	%2019-03-14 Created getGeneratedData function [by Jorrit Montijn]
	%2019-08-19 Added dblGaussNoise & boolDoublePeaked
	
	%% inputs
	if ~exist('dblGaussNoise','var') || isempty(dblGaussNoise)
		dblGaussNoise = 0;
	end
	if ~exist('boolDoublePeaked','var') || isempty(boolDoublePeaked)
		boolDoublePeaked = false;
	end
	
	%% generate preferred orientations
	vecPrefOri = rand(1,intPopN)*2*pi;
	intTrials = numel(vecTrialAngles);
	matResp = nan(intPopN,intTrials);
	if isscalar(vecKappa),vecKappa=repmat(vecKappa,[1 intPopN]);end
	for intN=1:intPopN
		%get mean tuning curve
		vecMeanR = circ_vmpdf(vecTrialAngles, vecPrefOri(intPopN), vecKappa(intPopN));
		if boolDoublePeaked
			vecMeanR = vecMeanR + circ_vmpdf(vecTrialAngles, vecPrefOri(intPopN)+deg2rad(180), vecKappa(intPopN));
		end
		
		%normalize to [10 dblDistDprime+10]
		vecMeanR = vecMeanR-min(vecMeanR);
		vecMeanR = dblDistDprime*(vecMeanR./max(vecMeanR));
		vecMeanR = vecMeanR+2;
		
		%generate responses by adding noise
		matResp(intN,:) = poissrnd(vecMeanR);
		
		%add Gaussian noise
		if dblGaussNoise > 0
		matResp(intN,:) = matResp(intN,:) + normrnd(0,dblGaussNoise,[1 intTrials]);
	end
end

