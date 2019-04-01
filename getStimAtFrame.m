function [vecOris,vecTrials] = getStimAtFrame(sIn,vecFrames)
	%getStimAtFrame Returns ori/trial at requested frame number
	%   Syntax: [vecOris,vecTrials] = getStimAtFrame(sIn,vecFrames)
	
	%check input type (ses or structStim)
	if isfield(sIn,'structStim')
		structStim = sIn.structStim;
	elseif isfield(sIn,'FrameOn')
		structStim = sIn;
	else
		error([mfilename ':StructureFormatNotDetected'],'Unknown structure format');
	end
	
	%pre-allocate
	vecOris = zeros(1,length(vecFrames));
	vecTrials = zeros(1,length(vecFrames));
	for intElement=1:length(vecFrames)
		%get orientation at frame
		lastStart = find(structStim.FrameOn < vecFrames(intElement),1,'last');
		lastStop = find(structStim.FrameOff < vecFrames(intElement),1,'last');
		if isempty(lastStart)
			vecOris(intElement) = 999;
			vecTrials(intElement) = 1;
		elseif isempty(lastStop)
			vecOris(intElement) = structStim.Orientation(1);
			vecTrials(intElement) = 1;
		else
			if lastStart > lastStop
				vecOris(intElement) = structStim.Orientation(lastStart);
			else
				vecOris(intElement) = 999;
			end
			vecTrials(intElement) = lastStart;
		end
	end
end

