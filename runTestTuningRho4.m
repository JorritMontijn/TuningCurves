%% set paramaters for data generation
clear all;close all;
strFigDir = 'D:\Data\ResultsOriMetric';
intN=100;
intNumPops = 50;

vecKappas = [20.342 9.106 5.174]; %FWHM: [15 22.5 30]
vecHzDiff = [0:0.5:3 4:20];
vecUniqueAngles = deg2rad(0:10:359);
intRandIters = 50;

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
%%
for intPlot=5:6
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
	dblFWHM = rad2deg(2*acos(1- [(1/dblKappa) * log(2)]))/2; %45 degs for kappa=9.106
	
	% pre-allocate
	matMeanDeltaPrimeBC = nan([numel(vecRep) numel(vecHzDiff) intNumPops]);
	matMeanTuningRho = nan([numel(vecRep) numel(vecHzDiff) intNumPops]);
	matMeanOSI = nan([numel(vecRep) numel(vecHzDiff) intNumPops]);
	matMeanOPI = nan([numel(vecRep) numel(vecHzDiff) intNumPops]);
	matCIDeltaPrimeBC = nan([numel(vecRep) numel(vecHzDiff) 2 intNumPops]);
	matCITuningRho = nan([numel(vecRep) numel(vecHzDiff) 2 intNumPops]);
	matCIOSI = nan([numel(vecRep) numel(vecHzDiff) 2 intNumPops]);
	matCIOPI = nan([numel(vecRep) numel(vecHzDiff) 2 intNumPops]);
	
	matMeanZscoreDeltaPrimeBC = nan([numel(vecRep) numel(vecHzDiff) intNumPops]);
	matMeanZscoreTuningRho = nan([numel(vecRep) numel(vecHzDiff) intNumPops]);
	matMeanZscoreOSI = nan([numel(vecRep) numel(vecHzDiff) intNumPops]);
	matMeanZscoreOPI = nan([numel(vecRep) numel(vecHzDiff) intNumPops]);
	matSdZscoreDeltaPrimeBC = nan([numel(vecRep) numel(vecHzDiff) intNumPops]);
	matSdZscoreTuningRho = nan([numel(vecRep) numel(vecHzDiff) intNumPops]);
	matSdZscoreOSI = nan([numel(vecRep) numel(vecHzDiff) intNumPops]);
	matSdZscoreOPI = nan([numel(vecRep) numel(vecHzDiff) intNumPops]);
	%%
	for intD=1:numel(vecHzDiff)
		dblHzDiff = vecHzDiff(intD);
		fprintf('Now at %.2f Hz (%d/%d) [%s]\n',dblHzDiff,intD,numel(vecHzDiff),getTime);
		for intC=1:numel(vecRep)
			%get nr of reps
			intRep = vecRep(intC);
			vecTrialAngles = repmat(vecUniqueAngles(:),[intRep 1])';
			
			% get number of pops
			for intPop=1:intNumPops
				% pre-allocate
				matShuffDeltaPrimeBC = nan(intN,intRandIters);
				matShuffTuningRho = nan(intN,intRandIters);
				matShuffOPI = nan(intN,intRandIters);
				matShuffOSI = nan(intN,intRandIters);
				
				%% get generated data
				vecKappa = dblKappa;%rand(1,intN)*1 + 10;
				[matResp,vecPrefOri] = getGeneratedData(intN,vecTrialAngles,vecKappa,dblHzDiff);
				
				% get bias-corrected delta prime
				vecDeltaPrimeBC = getDeltaPrime(matResp,vecTrialAngles,true);
				
				% get non-corrected delta prime
				vecTuningRho = getTuningRho(matResp,vecTrialAngles);
				
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
					matShuffTuningRho(:,intIter) = getTuningRho(matResp,vecShuffledTrialAngles,false);
					
					%get OSI
					matShuffOSI(:,intIter) = getOSI(matResp,vecShuffledTrialAngles);
					
					%get OPI
					matShuffOPI(:,intIter) = getOPI(matResp,vecShuffledTrialAngles);
				end
				
				%get means
				vecMeanShuffDeltaPrimeBC = nanmean(matShuffDeltaPrimeBC,2);
				vecSdShuffDeltaPrimeBC = nanstd(matShuffDeltaPrimeBC,[],2);
				vecMeanShuffTuningRho = nanmean(matShuffTuningRho,2);
				vecSdShuffTuningRho = nanstd(matShuffTuningRho,[],2);
				vecMeanShuffOSI= nanmean(matShuffOSI,2);
				vecSdShuffOSI= nanstd(matShuffOSI,[],2);
				vecMeanShuffOPI = nanmean(matShuffOPI,2);
				vecSdShuffOPI = nanstd(matShuffOPI,[],2);
				
				% get z-scores
				vecZscoreDeltaPrimeBC = (vecDeltaPrimeBC - vecMeanShuffDeltaPrimeBC) ./ vecSdShuffDeltaPrimeBC;
				vecZscoreTuningRho = (vecTuningRho - vecMeanShuffTuningRho) ./ vecSdShuffTuningRho;
				vecZscoreOSI = (vecOSI - vecMeanShuffOSI) ./ vecSdShuffOSI;
				vecZscoreOPI = (vecOPI - vecMeanShuffOPI) ./ vecSdShuffOPI;
				
				%get p-values
				vecPvalueDeltaPrimeBC = 2*normcdf(-abs(vecZscoreDeltaPrimeBC)) < 0.05;
				vecPvalueTuningRho = 2*normcdf(-abs(vecZscoreTuningRho)) < 0.05;
				vecPvalueOSI = 2*normcdf(-abs(vecZscoreOSI)) < 0.05;
				vecPvalueOPI = 2*normcdf(-abs(vecZscoreOPI)) < 0.05;
				
				%assign means + CI
				dblAlpha = 1 - 2*normcdf(-1); %one sd
				[phat,pci] = binofit(sum(vecPvalueDeltaPrimeBC),numel(vecPvalueDeltaPrimeBC),dblAlpha);
				matMeanDeltaPrimeBC(intC,intD,intPop) = phat;
				matCIDeltaPrimeBC(intC,intD,:,intPop) = pci;
				dblAlpha = 1 - 2*normcdf(-1); %one sd
				[phat,pci] = binofit(sum(vecPvalueTuningRho),numel(vecPvalueTuningRho),dblAlpha);
				matMeanTuningRho(intC,intD,intPop) = phat;
				matCITuningRho(intC,intD,:,intPop) = pci;
				[phat,pci] = binofit(sum(vecPvalueOSI),numel(vecPvalueOSI),dblAlpha);
				matMeanOSI(intC,intD,intPop) = phat;
				matCIOSI(intC,intD,:,intPop) = pci;
				[phat,pci] = binofit(sum(vecPvalueOPI),numel(vecPvalueOPI),dblAlpha);
				matMeanOPI(intC,intD,intPop) = phat;
				matCIOPI(intC,intD,:,intPop) = pci;
				
				%% z-scores
				%assign means
				matMeanZscoreDeltaPrimeBC(intC,intD,intPop) = nanmean(vecZscoreDeltaPrimeBC);
				matMeanZscoreTuningRho(intC,intD,intPop) = nanmean(vecZscoreTuningRho);
				matMeanZscoreOSI(intC,intD,intPop) = nanmean(vecZscoreOSI);
				matMeanZscoreOPI(intC,intD,intPop) = nanmean(vecZscoreOPI);
				
				%assign stds
				matSdZscoreDeltaPrimeBC(intC,intD,intPop) = nanstd(vecZscoreDeltaPrimeBC);
				matSdZscoreTuningRho(intC,intD,intPop) = nanstd(vecZscoreTuningRho);
				matSdZscoreOSI(intC,intD,intPop) = nanstd(vecZscoreOSI);
				matSdZscoreOPI(intC,intD,intPop) = nanstd(vecZscoreOPI);
			end
		end
	end
	
	%prep
	intSubPlot=intSubPlot+1;
	cellLegend = {'\delta''','\rho','OPI','OSI'};
	vecX = vecHzDiff;
	matX = repmat(vecX',[1 numel(cellLegend)]);
	
	%% plot z-score values
	matPlotY = cat(1,matMeanZscoreDeltaPrimeBC,matMeanZscoreTuningRho,matMeanZscoreOPI,matMeanZscoreOSI);
	matPlotY = mean(matPlotY,3)';
	vecEndVals = matPlotY(end,:);
	matPlotE = cat(1,matSdZscoreDeltaPrimeBC,matSdZscoreTuningRho,matSdZscoreOPI,matSdZscoreOSI);
	matPlotE = mean(matPlotE,3)';
	
	figure(h1);drawnow;
	subplot(2,3,intSubPlot)
	errorfill(matX,matPlotY,matPlotE);
	legend(cellLegend,'location','bestoutside')
	xlabel('d(Hz) between pref/non-pref')
	ylabel('Z-score (vs shuffled)')
	title(sprintf('%d rep,kappa=%.1f,FWHM=%.3f',vecRep,dblKappa,dblFWHM));
	fixfig
%	ylim([0 max(get(gca,'ylim'))]);
	
	%% plot p-value percentage
	if intSubPlot < 4
		intMaxHz = 4;
	else
		intMaxHz = 7;
	end
	figure(h2);drawnow;
	matX = matX(1:intMaxHz,:);
	matPlotY = cat(1,matMeanDeltaPrimeBC,matMeanTuningRho,matMeanOPI,matMeanOSI);
	matPlotY = mean(matPlotY,3)';
	matPlotY = matPlotY(1:intMaxHz,:);
	matPlotE = cat(1,matCIDeltaPrimeBC,matCITuningRho,matCIOPI,matCIOSI);
	matPlotE = mean(matPlotE,4);
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
strFig = sprintf('Bootstrapped_z-scores_%s',getDate);
export_fig(strcat(strFig,'.tif'));
print(gcf, '-dpdf', strcat(strFig,'.pdf'));
%% save 2
figure(h2);
drawnow;
strFig = sprintf('Bootstrapped_significance_%s',getDate);
export_fig(strcat(strFig,'.tif'));
print(gcf, '-dpdf', strcat(strFig,'.pdf'));

cd(strOldDir);
