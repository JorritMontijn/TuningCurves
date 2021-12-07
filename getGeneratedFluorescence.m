function [vecTimestamps,vecdFoF] = getGeneratedFluorescence(vecSpikeTimes,dblSamplingFreq,sIndicatorProps,boolQuick)
	%getGeneratedFluorescence Generates fluorescence data from spike times
	%   [vecTimestamps,vecdFoF] = getGeneratedFluorescence(vecSpikeTimes,dblSamplingFreq,sIndicatorProps,boolQuick)
	%
	%sIndicatorProps: structure with fields suppling properties of fluorescence indicator:
	%	- dblTimescale: scales indicator on/off half-lives (default: 1)
	%	- dblHalfT_On: rise half-life (default: 10/1000)
	%	- dblHalfT_Off decay half-life (default: 63/1000)
	%	- dblPeakMagnitude: peak magnitude in fraction over background (default: 0.05)
	%	- dblNoise: sd of Gaussian noise as fraction of background (default: 0.05)
	%
	%boolQuick: boolean, 
	%				true (default):
	%					uses filtering, assuming that all spike times occur at exactly the closest sample time
	%				false:
	%					uses subsample timing differences to calculate the indicator fluorescence
	%
	%Version History:
	%2021-11-29 Created getGeneratedFluorescence function [by Jorrit Montijn]
	%
	%GCaMP dynamics from Helassa et al., 2016 (sci reps): Design and mechanistic insight into ultrafast calcium indicators for monitoring intracellular calcium dynamics
	%t1/2 ON	t1/2 OFF
	%10+/-0.2	63?+/-6
	
	%% check inputs
	if dblSamplingFreq < 1
		warning([mfilename ':SampRateOrFreq'],'Sampling frequency is under 1 Hz, are you sure you are not supplying sampling interval?');
	end
	if ~exist('boolQuick','var') || isempty(boolQuick)
		boolQuick = true;
	end
	if ~exist('sIndicatorProps','var') || ~isstruct(sIndicatorProps)
		sIndicatorProps = struct;
	end
	dblTimescale = getOr(sIndicatorProps,'dblTimescale',3);
	dblHalfT_On = getOr(sIndicatorProps,'dblHalfT_On',dblTimescale*(10/1000));
	dblHalfT_Off = getOr(sIndicatorProps,'dblHalfT_Off',dblTimescale*(63/1000));
	dblPeakMagnitude = getOr(sIndicatorProps,'dblPeakMagnitude',5/100);
	dblNoise = getOr(sIndicatorProps,'dblNoise',5/100);
	
	%build sample times
	vecSampleEdges = 0:1/dblSamplingFreq:(max(vecSpikeTimes)+dblSamplingFreq);
	vecTimestamps = vecSampleEdges(2:end) - 0.5/dblSamplingFreq;
		
	%build kernel
	intKernelStart = ceil((dblHalfT_On*8)*dblSamplingFreq);
	intKernelEnd = ceil((dblHalfT_Off*8)*dblSamplingFreq);
	tau_on = dblHalfT_On/log(2);
	vecOn = dblPeakMagnitude*exp(-((1:intKernelStart)/dblSamplingFreq)/tau_on);
	tau_off = dblHalfT_Off/log(2);
	vecOff = dblPeakMagnitude*exp(-((1:intKernelEnd)/dblSamplingFreq)/tau_off);
	vecKernel = cat(2,vecOn(end:-1:1),dblPeakMagnitude,vecOff);
	intCenterT = intKernelStart + 1;
		
	%run
	if boolQuick
		%% ignore precise timings of spikes; assume they are exactly at the closest sample time
		%digitize spike times
		vecSpikeCounts = histcounts(vecSpikeTimes,vecSampleEdges);
		
		%apply
		vecTransientF = conv(vecSpikeCounts,vecKernel,'full');
		vecTransientF = vecTransientF(intCenterT:(numel(vecSpikeCounts)+intKernelStart));
	else
		%% take into account subsample location of spikes
		%create function
		fResp = @(vecSampT,dblSpikeT) double((vecSampT-dblSpikeT)<=0).*dblPeakMagnitude.*exp((vecSampT-dblSpikeT)/tau_on) + double((vecSampT-dblSpikeT)>0).*dblPeakMagnitude.*exp(-(vecSampT-dblSpikeT)/tau_off);
		intSampNum = numel(vecTimestamps);
		vecTransientF = zeros(size(vecTimestamps));
		for intSpike=1:numel(vecSpikeTimes)
			dblSpikeT = vecSpikeTimes(intSpike);
			
			intCenterSample = round(0.5+dblSpikeT*dblSamplingFreq);
			vecUseSamples = (intCenterSample-intKernelStart):(intCenterSample+intKernelEnd);
			vecUseSamples((vecUseSamples < 1) | (vecUseSamples > intSampNum)) = [];
			vecSampT = vecTimestamps(vecUseSamples);
			
			vecTransientF(vecUseSamples) = vecTransientF(vecUseSamples) + fResp(vecSampT,dblSpikeT);
		end
	end
	%% add noise & compute dF/F0
	%add noise
	vecF = 1 + vecTransientF + normrnd(0,dblNoise,size(vecTransientF));
	
	%transform to dF/F0
	[vecdFoF,vecF,vecF0] = calcLocaldFoF(vecF, dblSamplingFreq);
end
function [vecdFoF,vecF,vecF0] = calcLocaldFoF(vecF, dblSamplingFreq,boolSmooth)
	%set parameters
	dblFraction=0.5;
	dblSecWindowSize=30;
	if ~exist('boolSmooth','var') || isempty(boolSmooth)
		boolSmooth = false;
	end
	
	if boolSmooth
		%make smoothing kernel
		intKernelSteps = min([5 (round(dblSamplingFreq/4)*2-1)]); %kernel size has to be odd; set to be approximately half a second [or maximum of 5 to avoid over-smoothing]
		intKernelStepSize = 2/(intKernelSteps-1); %required size per step to get correct number of steps
		vecKernel = normpdf(-1:intKernelStepSize:1,0,1); %get gaussian kernel
		vecKernel = vecKernel / sum(vecKernel); %set integral to 1
		
		% smooth F
		vecF = imfilt(vecF(:)',vecKernel);
	end
	
	% calculate F0 window size
	dblWindowSecs = min( [ (0.4 * (length(vecF)/dblSamplingFreq)) dblSecWindowSize] ); %number of seconds for F0 baselining
	intWindowFrames = round(dblSamplingFreq*dblWindowSecs) ; %number of frames for F0 baselining
	
	boolLoop = (intWindowFrames*length(vecF)) > (1000)^3;
	if boolLoop
		% calculate F0 (baseline) per frame
		vecF0 = zeros(size(vecF)) ;
		sortF = sort( vecF( 1:intWindowFrames ) ) ;
		vecF0(1:intWindowFrames) = mean( sortF(1:round(intWindowFrames/2)) ) ;
		for i=(intWindowFrames+1):length(vecF)
			sortF = sort( vecF( i-intWindowFrames:i ) ) ;
			vecF0(i) = mean( sortF(1:round(intWindowFrames/2)) ) ;
		end

	else
		%calculate F0
		vecSelectBase=1:intWindowFrames; %base vector from which to build full selection matrix
		matSelect=repmat(vecSelectBase,[floor(intWindowFrames/2) 1]); %create static (non-shifting window) selection matrix for first part of trace
		intSizeSecond = (length(vecF)-intWindowFrames); %calculate size of second part of selection matrix
		matSelect=[matSelect;repmat(vecSelectBase,[intSizeSecond 1])+repmat((1:intSizeSecond)',[1 intWindowFrames])]; %add incrementally increasing selection trace matrix (with shifting window) to static first part
		matSelect=[matSelect;repmat(vecSelectBase,[ceil(intWindowFrames/2) 1])]; %create static (non-shifting window) selection matrix for last part of trace
		matSortedWindowTraces=sort(vecF(matSelect),2,'ascend');%select F values from trace with matSelect; then sort F values (ascending) per trace
		vecF0=mean(matSortedWindowTraces(:,1:round(intWindowFrames*dblFraction)),2)'; %take first half sorted values, so that lowest 50% are selected; then calculate mean per trace for each time point
	end
	%get dF/F
	vecdFoF = (vecF-vecF0)./vecF0; %calculate dF/F by subtracting F0 trace from F trace and dividing by F0 trace
end
