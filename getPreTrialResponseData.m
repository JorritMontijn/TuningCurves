function matPreTrialResponse = getPreTrialResponseData(ses,dblPreDur)
	%UNTITLED Summary of this function goes here
	%   Detailed explanation goes here
	
	%transform structure-based data to raw dFoF matrix
	intTrials = length(ses.structStim.FrameOn);
	intFrames = length(ses.neuron(1).dFoF);
	intNeurons = numel(ses.neuron);
	matActivity = zeros(intNeurons,intFrames);
	for intNeuron=1:intNeurons
		matActivity(intNeuron,:) = ses.neuron(intNeuron).dFoF;
	end
	
	%calculate mean per pre-trial
	matPreTrialResponse = zeros(intNeurons,intTrials);
	intDur = dblPreDur*round(ses.samplingFreq);
	for intTrial = 1:intTrials
		intFrameOn = ses.structStim.FrameOn(intTrial);
		intFrameStart = intFrameOn-intDur;
		matPreTrialResponse(:,intTrial) = mean(matActivity(:,intFrameStart:(intFrameOn-1)),2);
	end
end

