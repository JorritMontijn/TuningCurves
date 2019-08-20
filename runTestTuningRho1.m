%% set paramaters for data generation
strFigDir = 'D:\Data\ResultsOriMetric\';
intN=100;

vecKappa = [2.3 3 4 6 8 10 20 40 80 160 320];
vecUniqueAngles = deg2rad(0:5:359);
intRandIters = 1;

%new maximized figure
h1 = figure;
drawnow;
jFig = get(handle(gcf), 'JavaFrame');
jFig.setMaximized(true);
figure(gcf);
drawnow;

for intPlot=2
	% set params
	if intPlot == 1
		vecRep=10;
		dblHzDiff=5;%9.106;1;25
	elseif intPlot == 2
		vecRep=10;
		dblHzDiff=10;%9.106;1;25
	elseif intPlot == 3
		vecRep=100;
		dblHzDiff=5;%9.106;1;25
	elseif intPlot == 4
		vecRep=100;
		dblHzDiff=10;%9.106;1;25
		
	end
	
	% pre-allocate
	matMeanGammaShape = nan([numel(vecRep) numel(vecKappa)]);
	matSdGammaShape = nan([numel(vecRep) numel(vecKappa)]);
	
	
	for intD=1:numel(vecKappa)
		dblKappa = vecKappa(intD)
		for intC=1:numel(vecRep)
			%get nr of reps
			intRep = vecRep(intC);
			vecTrialAngles = repmat(vecUniqueAngles(:),[intRep 1])';
			vecTrialDegs = rad2deg(vecTrialAngles/2);
			
			% pre-allocate
			matTempGammaShape = nan(intN,intRandIters);
			for intIter=1:intRandIters
				%% get generated data
				[matResp,vecPrefOri] = getGeneratedData(intN,vecTrialAngles,dblKappa,dblHzDiff,0,1);
				
				%% get tuning parameter
				matTempGammaShape(:,intIter) = getTuningRho(matResp,vecTrialAngles);
				%[matTempGammaShape(:,intIter),sOut] = getVonMisesR2(matResp,vecTrialDegs);
				%vecFit = vonMisesSingleFitPX(sOut.matFittedParams(intNeuron,:),vecTrialDegs)
				%matResp
				
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
	vecFWHM = rad2deg(2*acos(1- [(1./(vecKappa)) * log(2)]))/2; %45 degs for kappa=9.106
	
	subplot(2,2,intPlot)
	errorfill(vecFWHM,matMeanGammaShape',matSdGammaShape');
	xlabel('FWHM (deg)')
	ylabel('Ori. selectivity')
	title(sprintf('%d repetitions, d(Hz)=%.1f',vecRep,dblHzDiff));
	fixfig
	%ylim([0 max(get(gca,'ylim'))]);
	
end
return
%% save 1
strOldDir = cd(strFigDir);
drawnow;
strFig = sprintf('TuningRho_FWHM_dependence');
export_fig(strcat(strFig,'.tif'));
print(gcf, '-dpdf', strcat(strFig,'.pdf'));
