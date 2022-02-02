function matResp = getRespMatSes(ses,vecNeurons,vecStims,structParams)
	%getRespMatSes Retrieves neuronal response for certain stimuli
	%	Syntax: matResp = getRespMatSes(ses,vecNeurons,vecStims,structParams)
	%   Input:
	%	- ses, session data
	%	- vecNeurons, vector of which neurons to include (optional)
	%	- vecStims, vector of which stimuli to include [-1 returns response
	%		outside stimulus presentations]; works well with cellSelect{}
	%		output vector (output from getSelectionVectors)  (optional)
	%	- structParams with fields: (optional)
	%	-	intStopOffset
	%	-	intStartOffset
	%	-	intPreBaselineRemoval
	%
	%	Output: 
	%	- matResp, 2D matrix containing neuronal response per stimulus per neuron:
	%		matResp(intNeuron,intStimPres)
	%
	%	Version history:
	%	1.0 - July 22 2013
	%	Created by Jorrit Montijn
	%	2.0 - May 19 2014
	%	Added support for preceding baseline subtraction [by JM]
	
	%check inputs
	if nargin < 4
		structParams = struct;
	end
	if nargin < 3 || isempty(vecStims)
		vecStims = 1:length(ses.structStim.FrameOff);
	end
	if nargin < 2 || isempty(vecNeurons)
		vecNeurons = 1:numel(ses.neuron);
	end
	
	%select frames
	if vecStims == -1
		%post-stimulus baseline
		vecStartFrames = ses.structStim.FrameOff;
		vecStopFrames = [ses.structStim.FrameOn(2:end) length(ses.neuron(1).dFoF)];
	elseif vecStims == -2
		%pre-stimulus baseline
		vecStartFrames = [round(ses.samplingFreq*3) ses.structStim.FrameOff(1:(end-1))];
		vecStopFrames = ses.structStim.FrameOn;
		vecStartFrames = round(max((vecStartFrames + vecStopFrames)/2,vecStartFrames-25));
	else
		%stimuli
		vecStartFrames = ses.structStim.FrameOn(vecStims);
		vecStopFrames = ses.structStim.FrameOff(vecStims);
	end
	
	%check if frame subset selection is requested
	if isfield(structParams,'intStopOffset')
		vecStopFrames = vecStartFrames + structParams.intStopOffset;
	end
	if isfield(structParams,'intStartOffset')
		vecStartFrames = vecStartFrames + structParams.intStartOffset;
	end
	if isfield(structParams,'intPreBaselineRemoval')
		intPreBaselineRemoval = structParams.intPreBaselineRemoval;%dblBaselineSecs
	else
		intPreBaselineRemoval = [];
	end
	
	%retrieve data
	if islogical(vecStims)
		intRepetitions = sum(vecStims);
	else
		intRepetitions = numel(vecStims);
	end
	%check if vector or scalar
	boolVec = length(vecNeurons) > 1;
	if boolVec
		if size(vecNeurons,2) == 1,vecNeurons = vecNeurons';end
		matResp = nan(max(vecNeurons),intRepetitions);
	else matResp = nan(1,intRepetitions);
	end
	
	%go through stims
	for intStimPres=1:length(vecStartFrames)
		intStartFrame = vecStartFrames(intStimPres);
		intStopFrame = vecStopFrames(intStimPres);
		if boolVec
			for intNeuron=vecNeurons
				if ~isempty(intPreBaselineRemoval),dblBaseline = mean(ses.neuron(intNeuron).dFoF((intStartFrame-intPreBaselineRemoval):(intStartFrame-1)));
				else dblBaseline = 0;end
				matResp(intNeuron,intStimPres) = mean(ses.neuron(intNeuron).dFoF(intStartFrame:intStopFrame))-dblBaseline;
			end; 
		else
			if ~isempty(intPreBaselineRemoval),dblBaseline = mean(ses.neuron(vecNeurons).dFoF((intStartFrame-intPreBaselineRemoval):(intStartFrame-1)));
			else dblBaseline = 0;end
			matResp(1,intStimPres) = mean(ses.neuron(vecNeurons).dFoF(intStartFrame:intStopFrame)) - dblBaseline;
		end
	end
end

