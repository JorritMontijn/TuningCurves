%% set paramaters for data generation
strFigDir = 'D:\Data\ResultsOriMetric\';
intN=100;
vecRep=3:100;
vecHzDiff = 5*[0 1 2 3 4];
vecUniqueAngles = deg2rad(0:45:359);
intRandIters = 1;
dblKappa=1;
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
vecX = vecRep(intStart:end);
matX = repmat(vecX',[1 5]);
matY_BC = matMeanDeltaPrimeBC(intStart:end,:);
matE_BC = matSdDeltaPrimeBC(intStart:end,:)*dblErrorbarSd;

cellLegend = vec2cell(vecHzDiff);
cellLegend = cellfun(@strcat,cellfill('d(Hz)=',size(cellLegend)),cellfun(@num2str,cellLegend,'UniformOutput',false),'UniformOutput',false);


subplot(2,2,1)
errorfill(matX,matY_BC,matE_BC);
legend(cellLegend,'location','bestoutside')
xlabel('# of repetitions')
ylabel('Ori. selectivity (\delta''_b_c)')
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

subplot(2,2,3)
errorfill(matX,matY_ICV,matE_ICV);
legend(cellLegend,'location','bestoutside')
xlabel('# of repetitions')
ylabel('Ori. selectivity (OPI)')
title('Ori. precision index (1 - Circ. var.)');
fixfig
ylim([0 max(get(gca,'ylim'))]);

%% OSI
matY_OI = matMeanOI(intStart:end,:);
matE_OI = matSdOI(intStart:end,:)*dblErrorbarSd;
subplot(2,2,4)
errorfill(matX,matY_OI,matE_OI);
legend(cellLegend,'location','bestoutside')
xlabel('# of repetitions')
ylabel('Ori. selectivity (OSI)','interpreter','none');
title('Ori. selectivity index ((P-NP)/P)')
fixfig
ylim([0 max(get(gca,'ylim'))]);



%% save
strOldDir = cd(strFigDir);
drawnow;
strFig = sprintf('Broadly_tuned_kappa=%d_FWHM=%.3f',dblKappa,dblFWHM);
export_fig(strcat(strFig,'.tif'));
print(gcf, '-dpdf', strcat(strFig,'.pdf'));
cd(strOldDir);

