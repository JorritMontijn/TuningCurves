function [vecRho_bc,vecRho,vecP] = getTuningRho(vecResp,vecAngles,boolBiasCorrection)
	%getTuningRho Calculates the tuning curve smoothness metric rho
	%	[vecRho_bc,vecRho,vecP] = getTuningRho(vecResp,vecAngles,boolBiasCorrection)
	%
	%Inputs:
	% - vecResp; [1 x Trial] response vector
	% - vecAngles; [1 x Trial] stimulus orientation vector  [in radians]
	% - boolBiasCorrection; [boolean] logical switch to use bias-correction
	%
	%Note: you can also supply a matrix of responses [Neuron x Trial]
	%instead of a vector, and the function will return a vector of rho
	%
	%Version History:
	%2019-04-12 Created tuning rho function [by Jorrit Montijn]
	
	
	%% check inputs
	if size(vecResp,2) ~= size(vecAngles,2)
		error([mfilename ':WrongInput'],'Response vector and angle vector are not the same length!');
	end
	if range(vecAngles) > 2*pi
		error([mfilename ':WrongInput'],'vecAngles is not in radians');
	end
	if ~exist('boolBiasCorrection','var') || isempty(boolBiasCorrection)
		boolBiasCorrection = true;
	end
	
	%% get range of continuum f & range of response r
	dblRangeF = range(vecAngles,2);
	dblRangeR = range(vecResp,2);
	
	%% calculate
	vecResp = bsxfun(@minus,vecResp,min(vecResp,[],2));
	[vecOriIdx,vecUniqueOris] = label2idx(vecAngles);
	intNumN = size(vecResp,1);
	intNumOri = max(vecOriIdx);
	intReps = numel(vecOriIdx)/intNumOri;
	intCombs = (intNumOri*(intNumOri-1))/2;
	matDeltaR = nan(intCombs,intNumN);
	matDeltaF = nan(intCombs,intNumN);
	intComb=0;
	for intOri1=1:intNumOri
		%get resp for ori 1
		vecRespOri1 = vecResp(:,vecOriIdx==intOri1);
		dblR1 = nanmean(vecRespOri1,2);
		dblF1 = vecUniqueOris(intOri1);
		for intOri2=(intOri1+1):intNumOri
			%get resp for ori 1
			vecRespOri2 = vecResp(:,vecOriIdx==intOri2);
			dblR2 = nanmean(vecRespOri2,2);
			dblF2 = vecUniqueOris(intOri2);
			dblDF = abs(circ_dist(dblF1,dblF2));
			dblDR = abs(dblR1 - dblR2);
			dblMuR = 1;%(dblR1+dblR2)/2;
			
			%calc lambda
			intComb=intComb+1;
			matDeltaR(intComb,:) = dblDR.*dblMuR;
			matDeltaF(intComb,:) = dblDF.*dblMuR;
		end
	end
	
	%% calc sd
	vecRho = nan(intNumN,1);
	vecRho_bc = nan(intNumN,1);
	vecP = nan(intNumN,1);
	for intN=1:intNumN
		%vecLambda2(intN)
		[bcR, p, T, df] = bcdistcorr(matDeltaR(:,intN),matDeltaF(:,intN));
		%bcR = corr(matDeltaR(:,intN),matDeltaF(:,intN));
		r2 = distcorr(matDeltaR(:,intN),matDeltaF(:,intN));
		vecRho(intN) = r2;
		vecRho_bc(intN) = bcR;
		vecP(intN) = p;
		%figure
		%scatter(matDeltaF(:,intN),matDeltaR(:,intN))
		%title(sprintf('r=%.3f; cr=%.3f, p=%.3f',bcR,r2,p))
	end
	
	if boolBiasCorrection
		%with bias-correction
		%dblDelta2Prime = dblDelta2Prime - (1 ./ (intReps.^(1.5)));
		%vecLambda = vecLambda  - log2(intNumOri)/(intReps.^(1.5));
			
	else
		%no bias-correction
	end
	
end

