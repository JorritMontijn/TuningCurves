function structOut = doVonMisesFit( ses, vecNeurons )
    %Fit a Von Mises distribution to the tuning curves of neurons to
    %determine the tuned orientation and the bandwidth (defined as the
    %angular distance between the tuned orientation and the half-width)
    
    if ~exist('vecNeurons','var') || isempty(vecNeurons)
		vecNeurons = 1:numel(ses.neuron);
    end
       
    %Get orientation angles
	vecStimTypeLookup = getOriListFromTrials(ses.structStim.Orientation);
    vecStimTypeLookup( vecStimTypeLookup == 999 ) = [];
    
    %Get tuning curves
    matTuning = getTuningCurve( ses, vecNeurons );
    sTuning = calcOriTuning( ses, vecNeurons );
    vecPrefAngle = circ_ang2rad( sTuning.vecPrefAngle(vecNeurons) );
    
    %Initialize variables
    vecTunedOri = nan(1,size(matTuning,1));
    vecAngleHalfWidth = nan(1,size(matTuning,1));
    
    for n = 1:size(matTuning,1)
        %Fit Von Mises function
        xvals = circ_ang2rad( vecStimTypeLookup ) * 2; %transform half circle into full circle angles
        yvals = squeeze(mean(matTuning(n,:,:),2));
        prefOri = vecPrefAngle(n);
        kappa = 1;
        
        %Set parameters for fit
        paramInit = [];
        paramInit(1) = prefOri*2;
        paramInit(2) = kappa;
        
        %Perform fit
        [paramBest, mse] = MLFit(@vonMises, paramInit, xvals, imnorm(yvals));
        YFitted = vonMises(circ_ang2rad(0:359), paramBest);
        
        %Get tuned orientation
        [dummy,vecTunedOri(n)] = max(YFitted);
        
        %Get angle at 1/sqrt(2) heigth
%         [~,vecAngleHalfWidth(n)] = min(abs( YFitted - (max(YFitted)*0.5) ));
        [dummy,vecAngleHalfWidth(n)] = min(abs( YFitted - (max(YFitted)*(1/sqrt(2))) ));
    end
    
    %Get angular distance between peak and half-width
    vecBandwidth = abs( circ_rad2ang( circ_dist(circ_ang2rad(vecTunedOri), circ_ang2rad(vecAngleHalfWidth)) ));
    
    %Transform circular angles back into half-circle angles
    structOut.vecTunedOri = vecTunedOri / 2;
    structOut.vecBandwidth = vecBandwidth / 2;
        

end

