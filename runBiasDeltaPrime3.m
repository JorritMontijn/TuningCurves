%% set paramaters for data generation
strFigDir = 'D:\Data\ResultsOriMetric\';
intN=100;
close all
vecHzDiff = [0:20];
vecUniqueAngles = deg2rad(0:45:359);
intRandIters = 1;

%new maximized figure
h1 = figure;
drawnow;
jFig = get(handle(gcf), 'JavaFrame');
jFig.setMaximized(true);
figure(gcf);
drawnow;
%new maximized figure
h2 = figure;
drawnow;
jFig = get(handle(gcf), 'JavaFrame');
jFig.setMaximized(true);
figure(gcf);
drawnow;
%set plot 0
intSubPlot=0;

for intPlot=1:4
	% set params
	if intPlot == 1
		vecRep=100;
		dblKappa=25;%9.106;1;25
	elseif intPlot == 2
		vecRep=100;
		dblKappa=1;%9.106;1;25
	elseif intPlot == 3
		vecRep=20;
		dblKappa=25;%9.106;1;25
	elseif intPlot == 4
		vecRep=20;
		dblKappa=1;%9.106;1;25
		
	end
	dblFWHM = rad2deg(2*acos(1- [(1/dblKappa) * log(2)])); %45 degs for kappa=9.106
	
	% pre-allocate
	matMeanDeltaPrimeBC = nan([numel(vecRep) numel(vecHzDiff)]);
	matMeanDeltaPrime = nan([numel(vecRep) numel(vecHzDiff)]);
	matMeanOI = nan([numel(vecRep) numel(vecHzDiff)]);
	matMeanICV = nan([numel(vecRep) numel(vecHzDiff)]);
	matSdDeltaPrimeBC = nan([numel(vecRep) numel(vecHzDiff)]);
	matSdDeltaPrime = nan([numel(vecRep) numel(vecHzDiff)]);
	matSdOI = nan([numel(vecRep) numel(vecHzDiff)]);
	matSdICV = nan([numel(vecRep) numel(vecHzDiff)]);
	
	
	for intD=1:numel(vecHzDiff)
		dblHzDiff = vecHzDiff(intD)
		for intC=1:numel(vecRep)
			%get nr of reps
			intRep = vecRep(intC);
			vecTrialAngles = repmat(vecUniqueAngles(:),[intRep 1])';
			
			% pre-allocate
			matTempDeltaPrimeBC = nan(intN,intRandIters);
			matTempDeltaPrime = nan(intN,intRandIters);
			matTempICV = nan(intN,intRandIters);
			matTempOI = nan(intN,intRandIters);
			for intIter=1:intRandIters
				%% get generated data
				vecKappa = dblKappa;%rand(1,intN)*1 + 10;
				[matResp,vecPrefOri] = getGeneratedData(intN,vecTrialAngles,vecKappa,dblHzDiff);
				
				%% get bias-corrected delta prime
				matTempDeltaPrimeBC(:,intIter) = getDeltaPrime(matResp,vecTrialAngles,true);
				
				%% get non-corrected delta prime
				matTempDeltaPrime(:,intIter) = getDeltaPrime(matResp,vecTrialAngles,false);
				
				
				%OI = (Rpref_ori ? Rorth)/Rpref_ori
				vecAngleIdx = label2idx(vecTrialAngles);
				vecPref = nan(1,intN);
				vecNonPref = nan(1,intN);
				for intN=1:size(matResp,1)
					vecR = accumarray(vecAngleIdx',matResp(intN,:));
					%% get OI
					[dblPref,intIdx] = max(vecR);
					dblPrefAngle = vecUniqueAngles(intIdx);
					[dblOrthAngle,intOrthIdx] = min(abs(circ_dist(vecUniqueAngles,dblPrefAngle+deg2rad(90))));
					if intOrthIdx == 0,intOrthIdx=numel(vecUniqueAngles);end
					matTempOI(intN,intIter) = (dblPref - vecR(intOrthIdx)) ./ dblPref;
					
					%% get 1 - circ_var
					matTempICV(intN,intIter) = 1 - circ_var(vecUniqueAngles',vecR);
				end
				
				
			end
			
			%get means
			vecMeanDeltaPrimeBC = nanmean(matTempDeltaPrimeBC,2);
			vecMeanDeltaPrime = nanmean(matTempDeltaPrime,2);
			vecMeanOI= nanmean(matTempOI,2);
			vecMeanICV = nanmean(matTempICV,2);
			
			%assign means
			matMeanDeltaPrimeBC(intC,intD) = nanmean(vecMeanDeltaPrimeBC);
			matMeanDeltaPrime(intC,intD) = nanmean(vecMeanDeltaPrime);
			matMeanOI(intC,intD) = nanmean(vecMeanOI);
			matMeanICV(intC,intD) = nanmean(vecMeanICV);
			
			%assign stds
			matSdDeltaPrimeBC(intC,intD) = nanstd(vecMeanDeltaPrimeBC);
			matSdDeltaPrime(intC,intD) = nanstd(vecMeanDeltaPrime);
			matSdOI(intC,intD) = nanstd(vecMeanOI);
			matSdICV(intC,intD) = nanstd(vecMeanICV);
		end
	end
	
	
	%% plot normal
	figure(h1);drawnow;
	intSubPlot=intSubPlot+1;
	dblErrorbarSd = 1;
	intStart = 1;
	intD = 1;
	vecX = vecHzDiff;
	matX = repmat(vecX',[1 3]);
	matPlotY = cat(1,matMeanDeltaPrime,matMeanICV,matMeanOI)';
	matPlotE = cat(1,matSdDeltaPrime,matSdICV,matSdOI)';
	cellLegend = {'\delta''','OPI','OSI'};
	
	subplot(2,2,intSubPlot)
	errorfill(matX,matPlotY,matPlotE);
	legend(cellLegend,'location','bestoutside')
	xlabel('d(Hz) between pref/non-pref')
	ylabel('Ori. selectivity')
	title(sprintf('%d repetitions, kappa=%.1f, FWHM=%.3f',vecRep,dblKappa,dblFWHM));
	fixfig
	ylim([0 max(get(gca,'ylim'))]);
	
	%% plot normalized
	matPlotY = cat(1,matMeanDeltaPrime,matMeanICV,matMeanOI)';
	vecEndVals = matPlotY(end,:);
	matPlotY = bsxfun(@rdivide,matPlotY,vecEndVals);
	matPlotE = cat(1,matSdDeltaPrime,matSdICV,matSdOI)';
	matPlotE = bsxfun(@rdivide,matPlotE,vecEndVals);
	
	figure(h2);drawnow;
	subplot(2,2,intSubPlot)
	errorfill(matX,matPlotY,matPlotE);
	legend(cellLegend,'location','bestoutside')
	xlabel('d(Hz) between pref/non-pref')
	ylabel('Normalized selectivity')
	title(sprintf('%d repetitions, kappa=%.1f, FWHM=%.3f',vecRep,dblKappa,dblFWHM));
	fixfig
	ylim([0 max(get(gca,'ylim'))]);
end

%% save 1
strOldDir = cd(strFigDir);
figure(h1);
drawnow;
strFig = sprintf('dMu_dependence');
export_fig(strcat(strFig,'.tif'));
print(gcf, '-dpdf', strcat(strFig,'.pdf'));
%% save 2
figure(h2);
drawnow;
strFig = sprintf('Norm_dMu_dependence');
export_fig(strcat(strFig,'.tif'));
print(gcf, '-dpdf', strcat(strFig,'.pdf'));

cd(strOldDir);
