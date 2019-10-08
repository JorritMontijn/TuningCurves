%% set paramaters for data generation
clear all;
strFigDir = 'D:\Data\ResultsOriMetric\';
intMaxN=100;
intReps=20;
dblDistDprime = 0;%5*[0 1 2 3 4];
vecUniqueAngles = deg2rad(0:45:359);
intRandIters = 1;
dblKappa=25;%9.106;1;25
dblFWHM = rad2deg(2*acos(1- [(1/dblKappa) * log(2)])); %45 degs for kappa=9.106

%get nr of reps
vecTrialAngles = repmat(vecUniqueAngles(:),[intReps 1])';
vecTrialIdx = label2idx(vecTrialAngles);
intOris = numel(vecUniqueAngles);

%% get generated data
[matResp,vecPrefOri] = getGeneratedData(intMaxN,vecTrialAngles,dblKappa,dblDistDprime,0.5);
matResp = bsxfun(@minus,matResp,min(matResp,[],2));
	
%% get non-corrected delta prime
vecDeltaPrime = getDeltaPrime(matResp,vecTrialAngles,false);
[vecOPI,vecROPI] = getOPI(matResp,vecTrialAngles);
[vecOSI,vecROSI,intIdx,intOrthIdx] = getOSI(matResp,vecTrialAngles);
matR = nan(intOris,intReps,intMaxN);

%% transform
for intN=1:intMaxN
	for intOri=1:intOris
		matR(intOri,:,intN) = matResp(intN,vecTrialIdx==intOri);
	end
end

%% plot neuron
[d,intPlotN] = max(vecOPI);
dblDP = vecDeltaPrime(intPlotN);
dblOSI = vecOSI(intPlotN);
dblOPI = vecOPI(intPlotN);

vecOriDeg=rad2deg(vecUniqueAngles);
matPlotR = matR(:,:,intPlotN);
vecMean = mean(matPlotR,2);
vecSd = std(matPlotR,[],2);
[dblOSI2,vecROSI2,intIdx2,intOrthIdx2] = getOSI(vecMean',vecUniqueAngles);
[dblOPI2,vecROPI2] = getOPI(vecMean',vecUniqueAngles);

subplot(2,2,1)
plot(vecOriDeg,vecMean);
xlabel('Stim Ori (deg)')
ylabel('Mean spiking rate (Hz)')
title(sprintf('OSI=%.3f;OPI=%.3f',dblOSI2,dblOPI2))
fixfig;

subplot(2,2,2)
scatter(flat(repmat(vecOriDeg',[1 intReps])),flat(matPlotR+min(matPlotR(:))))
xlabel('Stim Ori (deg)')
ylabel('Spiking rate per trial (Hz)')
title(sprintf('OSI=%.3f;OPI=%.3f;d''=%.3f',dblOSI,dblOPI,dblDP))
fixfig;

%%
return
%new maximized figure
figure
drawnow;
jFig = get(handle(gcf), 'JavaFrame');
jFig.setMaximized(true);
figure(gcf);
drawnow;

%% delta' bc
dblErrorbarSd = 1;
intStart = 1;
intD = 1;
vecX = intReps(intStart:end);
matX = repmat(vecX',[1 5]);
matY_BC = matMeanDeltaPrimeBC(intStart:end,:);
matE_BC = matSdDeltaPrimeBC(intStart:end,:)*dblErrorbarSd;

cellLegend = vec2cell(dblDistDprime);
cellLegend = cellfun(@strcat,cellfill('d''=',size(cellLegend)),cellfun(@num2str,cellLegend,'UniformOutput',false),'UniformOutput',false);


subplot(2,2,1)
errorfill(matX,matY_BC,matE_BC);
legend(cellLegend,'location','bestoutside')
xlabel('# of repetitions')
ylabel('Ori. selectivity (\delta'')')
title('Bias-corrected orientation metric \delta''');
fixfig
ylim([0 max(get(gca,'ylim'))]);

%% delta' non-bc
matY = matMeanDeltaPrime(intStart:end,:);
matE = matSdDeltaPrime(intStart:end,:)*dblErrorbarSd;
subplot(2,2,2)
errorfill(matX,matY,matE);
legend(cellLegend,'location','bestoutside')
xlabel('# of repetitions')
ylabel('Ori. selectivity (\delta'')')
title('Non-bias-corrected orientation metric \delta''');
fixfig
ylim([0 max(get(gca,'ylim'))]);

%% ICV
matY_ICV = matMeanICV(intStart:end,:);
matE_ICV = matSdICV(intStart:end,:)*dblErrorbarSd;

cellLegend = vec2cell(dblDistDprime);
cellLegend = cellfun(@strcat,cellfill('d''=',size(cellLegend)),cellfun(@num2str,cellLegend,'UniformOutput',false),'UniformOutput',false);
subplot(2,2,3)
errorfill(matX,matY_ICV,matE_ICV);
legend(cellLegend,'location','bestoutside')
xlabel('# of repetitions')
ylabel('Ori. selectivity (1 - CV)')
title('Inverse of circular variance');
fixfig
ylim([0 max(get(gca,'ylim'))]);

%% OI
matY_OI = matMeanOI(intStart:end,:);
matE_OI = matSdOI(intStart:end,:)*dblErrorbarSd;
subplot(2,2,4)
errorfill(matX,matY_OI,matE_OI);
legend(cellLegend,'location','bestoutside')
xlabel('# of repetitions')
ylabel('Ori. selectivity ((P-NP)/P)')
title('Orientation index, (R_pref - R_nonpref)/R_pref','interpreter','none');
fixfig
ylim([0 max(get(gca,'ylim'))]);


