%% set paramaters for data generation
strFigDir = 'D:\Data\ResultsOriMetric\';
intN=100;

vecHzDiff = [0:20];
vecUniqueAngles = deg2rad(0:20:359);
intRandIters = 1;

%new maximized figure
h1 = figure;
drawnow;
jFig = get(handle(gcf), 'JavaFrame');
jFig.setMaximized(true);
figure(gcf);
drawnow;

for intPlot=1:4
	% set params
	if intPlot == 1
		vecRep=100;
		dblKappa=10;%9.106;1;25
	elseif intPlot == 2
		vecRep=100;
		dblKappa=1;%9.106;1;25
	elseif intPlot == 3
		vecRep=20;
		dblKappa=10;%9.106;1;25
	elseif intPlot == 4
		vecRep=20;
		dblKappa=1;%9.106;1;25
		
	end
	dblFWHM = rad2deg(2*acos(1- [(1/dblKappa) * log(2)]))/2; %45 degs for kappa=9.106
	
	% pre-allocate
	matMeanGammaShape = nan([numel(vecRep) numel(vecHzDiff)]);
	matSdGammaShape = nan([numel(vecRep) numel(vecHzDiff)]);
	
	
	for intD=1:numel(vecHzDiff)
		dblHzDiff = vecHzDiff(intD)
		for intC=1:numel(vecRep)
			%get nr of reps
			intRep = vecRep(intC);
			vecTrialAngles = repmat(vecUniqueAngles(:),[intRep 1])';
			vecTrialDegs = rad2deg(vecTrialAngles/2);
			
			% pre-allocate
			matTempGammaShape = nan(intN,intRandIters);
			for intIter=1:intRandIters
				%% get generated data
				vecKappa = dblKappa;%rand(1,intN)*1 + 10;
				[matResp,vecPrefOri] = getGeneratedData(intN,vecTrialAngles,vecKappa,dblHzDiff);
				
				%% get bias-corrected delta prime
				matTempGammaShape(:,intIter) = getTuningRho(matResp,vecTrialAngles);
				%[matTempGammaShape(:,intIter),sOut] = getVonMisesR2(matResp,vecTrialDegs);
			
				
			end
			
			%get means
			vecMeanGammaShape = nanmean(matTempGammaShape,2);
			
			%assign means
			matMeanGammaShape(intC,intD) = nanmean(vecMeanGammaShape);
			
			%assign stds
			matSdGammaShape(intC,intD) = nanstd(vecMeanGammaShape);
			
		end
	end
	
	
	%% plot normal
	figure(h1);drawnow;
	dblErrorbarSd = 1;
	intStart = 1;
	intD = 1;
	
	subplot(2,2,intPlot)
	errorfill(vecHzDiff,matMeanGammaShape',matSdGammaShape');
	xlabel('d(Hz) between pref/non-pref')
	ylabel('Ori. selectivity')
	title(sprintf('%d repetitions, kappa=%.1f, FWHM=%.3f',vecRep,dblKappa,dblFWHM));
	fixfig
	%ylim([0 max(get(gca,'ylim'))]);
	
end
return
%% save 1
strOldDir = cd(strFigDir);
%figure(h1);
drawnow;
strFig = sprintf('TuningRho_dMu_dependence');
export_fig(strcat(strFig,'.tif'));
print(gcf, '-dpdf', strcat(strFig,'.pdf'));
