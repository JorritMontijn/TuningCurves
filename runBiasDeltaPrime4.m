%% set paramaters for data generation
clear all;close all;
strFigDir = 'D:\Data\ResultsOriMetricDeltaPrime\';
intN=100;

vecHzDiff = [0:0.5:3 4:20];
vecUniqueAngles = deg2rad(0:45:359);
intRandIters = 100;

%new maximized figure
h1 = figure;
drawnow;
jFig = get(handle(gcf), 'JavaFrame');
jFig.setMaximized(true);
figure(gcf);
drawnow;
h2 = figure;
drawnow;
jFig = get(handle(gcf), 'JavaFrame');
jFig.setMaximized(true);
figure(gcf);
drawnow;

%set plot 0
intSubPlot=0;

for intPlot=1:6
	% set params
	if intPlot == 1
		vecRep=100;
		dblKappa=25;%9.106;1;25
	elseif intPlot == 2
		vecRep=100;
		dblKappa=9.106;
	elseif intPlot == 3
		vecRep=100;
		dblKappa=1;%9.106;1;25
	elseif intPlot == 4
		vecRep=20;
		dblKappa=25;%9.106;1;25
	elseif intPlot == 5
		vecRep=20;
		dblKappa=9.106;%9.106;1;25
	elseif intPlot == 6
		vecRep=20;
		dblKappa=1;%9.106;1;25
		
	end
	dblFWHM = rad2deg(2*acos(1- [(1/dblKappa) * log(2)])); %45 degs for kappa=9.106
	
	% pre-allocate
	matMeanDeltaPrimeBC = nan([numel(vecRep) numel(vecHzDiff)]);
	matMeanDeltaPrime = nan([numel(vecRep) numel(vecHzDiff)]);
	matMeanOSI = nan([numel(vecRep) numel(vecHzDiff)]);
	matMeanOPI = nan([numel(vecRep) numel(vecHzDiff)]);
	matCIDeltaPrimeBC = nan([numel(vecRep) numel(vecHzDiff) 2]);
	matCIDeltaPrime = nan([numel(vecRep) numel(vecHzDiff) 2]);
	matCIOSI = nan([numel(vecRep) numel(vecHzDiff) 2]);
	matCIOPI = nan([numel(vecRep) numel(vecHzDiff) 2]);
	
	matMeanZscoreDeltaPrimeBC = nan([numel(vecRep) numel(vecHzDiff)]);
	matMeanZscoreDeltaPrime = nan([numel(vecRep) numel(vecHzDiff)]);
	matMeanZscoreOSI = nan([numel(vecRep) numel(vecHzDiff)]);
	matMeanZscoreOPI = nan([numel(vecRep) numel(vecHzDiff)]);
	matSdZscoreDeltaPrimeBC = nan([numel(vecRep) numel(vecHzDiff)]);
	matSdZscoreDeltaPrime = nan([numel(vecRep) numel(vecHzDiff)]);
	matSdZscoreOSI = nan([numel(vecRep) numel(vecHzDiff)]);
	matSdZscoreOPI = nan([numel(vecRep) numel(vecHzDiff)]);
	
	for intD=1:numel(vecHzDiff)
		dblHzDiff = vecHzDiff(intD)
		for intC=1:numel(vecRep)
			%get nr of reps
			intRep = vecRep(intC);
			vecTrialAngles = repmat(vecUniqueAngles(:),[intRep 1])';
			
			% pre-allocate
			matShuffDeltaPrimeBC = nan(intN,intRandIters);
			matShuffDeltaPrime = nan(intN,intRandIters);
			matShuffOPI = nan(intN,intRandIters);
			matShuffOSI = nan(intN,intRandIters);
			
			%% get generated data
			vecKappa = dblKappa;%rand(1,intN)*1 + 10;
			[matResp,vecPrefOri] = getGeneratedData(intN,vecTrialAngles,vecKappa,dblHzDiff);
			
			% get bias-corrected delta prime
			vecDeltaPrimeBC = getDeltaPrime(matResp,vecTrialAngles,true);
			
			% get non-corrected delta prime
			vecDeltaPrime = getDeltaPrime(matResp,vecTrialAngles,false);
			
			%get OSI
			vecOSI = getOSI(matResp,vecTrialAngles);
			
			%get OPI
			vecOPI = getOPI(matResp,vecTrialAngles);
			
			
			%% get shuffled bootstrap
			for intIter=1:intRandIters
				vecShuffledTrialAngles = vecTrialAngles(randperm(numel(vecTrialAngles)));
				% get bias-corrected delta prime
				matShuffDeltaPrimeBC(:,intIter) = getDeltaPrime(matResp,vecShuffledTrialAngles,true);
				
				% get non-corrected delta prime
				matShuffDeltaPrime(:,intIter) = getDeltaPrime(matResp,vecShuffledTrialAngles,false);
				
				%get OSI
				matShuffOSI(:,intIter) = getOSI(matResp,vecShuffledTrialAngles);
				
				%get OPI
				matShuffOPI(:,intIter) = getOPI(matResp,vecShuffledTrialAngles);
			end
			
			%get means
			vecMeanShuffDeltaPrimeBC = nanmean(matShuffDeltaPrimeBC,2);
			vecSdShuffDeltaPrimeBC = nanstd(matShuffDeltaPrimeBC,[],2);
			vecMeanShuffDeltaPrime = nanmean(matShuffDeltaPrime,2);
			vecSdShuffDeltaPrime = nanstd(matShuffDeltaPrime,[],2);
			vecMeanShuffOSI= nanmean(matShuffOSI,2);
			vecSdShuffOSI= nanstd(matShuffOSI,[],2);
			vecMeanShuffOPI = nanmean(matShuffOPI,2);
			vecSdShuffOPI = nanstd(matShuffOPI,[],2);
			
			% get z-scores
			vecZscoreDeltaPrimeBC = (vecDeltaPrimeBC - vecMeanShuffDeltaPrimeBC) ./ vecSdShuffDeltaPrimeBC;
			vecZscoreDeltaPrime = (vecDeltaPrime - vecMeanShuffDeltaPrime) ./ vecSdShuffDeltaPrime;
			vecZscoreOSI = (vecOSI - vecMeanShuffOSI) ./ vecSdShuffOSI;
			vecZscoreOPI = (vecOPI - vecMeanShuffOPI) ./ vecSdShuffOPI;
			
			%get p-values
			vecPvalueDeltaPrime = 2*normcdf(-abs(vecZscoreDeltaPrime)) < 0.05;
			vecPvalueOSI = 2*normcdf(-abs(vecZscoreOSI)) < 0.05;
			vecPvalueOPI = 2*normcdf(-abs(vecZscoreOPI)) < 0.05;
			
			%assign means + CI
			dblAlpha = 1 - 2*normcdf(-1); %one sd
			[phat,pci] = binofit(sum(vecPvalueDeltaPrime),numel(vecPvalueDeltaPrime),dblAlpha);
			matMeanDeltaPrimeBC(intC,intD) = phat;
			matCIDeltaPrimeBC(intC,intD,:) = pci;
			[phat,pci] = binofit(sum(vecPvalueOSI),numel(vecPvalueOSI),dblAlpha);
			matMeanOSI(intC,intD) = phat;
			matCIOSI(intC,intD,:) = pci;
			[phat,pci] = binofit(sum(vecPvalueOPI),numel(vecPvalueOPI),dblAlpha);
			matMeanOPI(intC,intD) = phat;
			matCIOPI(intC,intD,:) = pci;
			
			%% z-scores
			%assign means
			matMeanZscoreDeltaPrimeBC(intC,intD) = nanmean(vecZscoreDeltaPrime);
			matMeanZscoreOSI(intC,intD) = nanmean(vecZscoreOSI);
			matMeanZscoreOPI(intC,intD) = nanmean(vecZscoreOPI);
			
			%assign stds
			matSdZscoreDeltaPrimeBC(intC,intD) = nanstd(vecZscoreDeltaPrime);
			matSdZscoreOSI(intC,intD) = nanstd(vecZscoreOSI);
			matSdZscoreOPI(intC,intD) = nanstd(vecZscoreOPI);
		end
	end
	
	%prep
	intSubPlot=intSubPlot+1;
	vecX = vecHzDiff;
	matX = repmat(vecX',[1 3]);
	cellLegend = {'\delta''','OPI','OSI'};
	%% plot z-score values
	matPlotY = cat(1,matMeanZscoreDeltaPrimeBC,matMeanZscoreOPI,matMeanZscoreOSI)';
	vecEndVals = matPlotY(end,:);
	matPlotE = cat(1,matSdZscoreDeltaPrimeBC,matSdZscoreOPI,matSdZscoreOSI)';
	
	figure(h1);drawnow;
	subplot(2,3,intSubPlot)
	errorfill(matX,matPlotY,matPlotE);
	legend(cellLegend,'location','bestoutside')
	xlabel('d(Hz) between pref/non-pref')
	ylabel('Z-score (vs shuffled)')
	title(sprintf('%d rep,kappa=%.1f,FWHM=%.3f',vecRep,dblKappa,dblFWHM));
	fixfig
	ylim([0 max(get(gca,'ylim'))]);
	
	%% plot p-value percentage
	if intSubPlot < 4
		intMaxHz = 4;
	else
		intMaxHz = 7;
	end
	figure(h2);drawnow;
	matX = matX(1:intMaxHz,:);
	matPlotY = cat(1,matMeanDeltaPrimeBC,matMeanOPI,matMeanOSI)';
	matPlotY = matPlotY(1:intMaxHz,:);
	matPlotE = cat(1,matCIDeltaPrimeBC,matCIOPI,matCIOSI);
	matPlotE = matPlotE(:,1:intMaxHz,:);
	matPlotE = bsxfun(@minus,matPlotY,permute(matPlotE,[2 1 3]));
	matPlotE(:,:,[1 2]) = abs(matPlotE(:,:,[2 1]));
	
	subplot(2,3,intSubPlot)
	errorfill(matX,matPlotY,matPlotE);
	legend(cellLegend,'location','bestoutside')
	xlabel('d(Hz) between pref/non-pref')
	ylabel('Fraction of sign. tuned neurons')
	title(sprintf('%d rep,kappa=%.1f,FWHM=%.3f',vecRep,dblKappa,dblFWHM));
	fixfig
	ylim([0 max(get(gca,'ylim'))]);
	
end

%% save 1
strOldDir = cd(strFigDir);
figure(h1);
drawnow;
strFig = sprintf('Bootstrapped_z-scores');
export_fig(strcat(strFig,'.tif'));
print(gcf, '-dpdf', strcat(strFig,'.pdf'));
%% save 2
figure(h2);
drawnow;
strFig = sprintf('Bootstrapped_significance');
export_fig(strcat(strFig,'.tif'));
print(gcf, '-dpdf', strcat(strFig,'.pdf'));

cd(strOldDir);
