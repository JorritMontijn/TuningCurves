%% set paramaters for data generation
intN=10;
vecRep=3:100;
vecDistDprime = 5*[0 1 2 3 4];
vecUniqueAngles = deg2rad(0:20:359);
intRandIters = 1;
dblKappa=4;%9.106;1;25
dblFWHM = rad2deg(2*acos(1- [(1/dblKappa) * log(2)]))/2; %45 degs for kappa=9.106

% pre-allocate
matMeanGammaShape = nan([numel(vecRep) numel(vecDistDprime)]);
matSdGammaShape = nan([numel(vecRep) numel(vecDistDprime)]);


for intD=1:numel(vecDistDprime)
	dblDistDprime = vecDistDprime(intD)
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
			[matResp,vecPrefOri] = getGeneratedData(intN,vecTrialAngles,vecKappa,dblDistDprime);
			
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
matY_BC = matMeanGammaShape(intStart:end,:);
matE_BC = matSdGammaShape(intStart:end,:)*dblErrorbarSd;

cellLegend = vec2cell(vecDistDprime);
cellLegend = cellfun(@strcat,cellfill('d''=',size(cellLegend)),cellfun(@num2str,cellLegend,'UniformOutput',false),'UniformOutput',false);


errorfill(matX,matY_BC,matE_BC);
legend(cellLegend,'location','bestoutside')
xlabel('# of repetitions')
ylabel('Ori. selectivity (\delta'')')
title('Tuning metric rho');
fixfig
ylim([-0.1 max(get(gca,'ylim'))]);
