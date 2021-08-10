function [matTE,vecWindowBinCenters] = getRespMat(vecTime,vecVals,vecEvents,vecWindow)
	%getRespMat Retrieves [time x event] matrix
	%	Syntax: matTE = getRespMat(vecTime,vecVals,vecEvents,vecWindow)
	%   Input:
	%	- vecTime, vector of time stamps
	%	- vecVals, vector of values
	%	- vecEvents, vector of event times
	%	- vecWindow, 2-element vector specifying start and stop times;
	%					or vector with bin edges
	%	Output:
	%	- matTE, 2D matrix containing binned value per time point per event:
	%		matTE(intTimeBin,intEvent)
	%	- vecWindowBinCenters, vector containing timestamps per bin
	%
	%	Version history:
	%	1.0 - 2021 August 10
	%	Created by Jorrit Montijn
	
	%% input is trace
	%get window
	if ~exist('vecWindow','var') || isempty(vecWindow)
		vecWindow = [-1 3];
	end
	
	%get event times
	intEvents = numel(vecEvents);
	vecEventStarts = min(vecWindow) + vecEvents;
	vecEventStops = max(vecWindow) + vecEvents;
	
	if numel(vecWindow) == 2
		%get window variables
		intWindowSize = 1 + find(vecTime >= vecEventStops(2),1) - find(vecTime >= vecEventStarts(2),1);
		vecWindowBinCenters = (0:(intWindowSize-1))/(intWindowSize-1);
		vecWindowBinCenters = (vecWindowBinCenters * range(vecWindow)) + vecWindow(1);
	else
		vecWindowBinCenters = vecWindow(1:(end-1)) + diff(vecWindow)/2;
	end
	
	%use simple trial loop
	intWindowSize = numel(vecWindowBinCenters);
	matTE = nan(intEvents,intWindowSize);
	for intEvent=1:intEvents
		if numel(vecWindow) == 2
			%retrieve target entries
			vecAssignPoints = 1:intWindowSize;
			intStart = find(vecTime >= vecEventStarts(intEvent),1);
			intStop = find(vecTime >= vecEventStops(intEvent),1);
			if isempty(intStop) %out-of-bounds at end
				intStop=numel(vecTime);
				vecUsePoints = intStart:intStop;
				vecAssignPoints((numel(vecUsePoints)+1):end) = []; %remove out-of-bounds entries
			end
			if intStart == 1 %out-of-bounds at start
				vecUsePoints = intStart:intStop;
				intDeleteUpTo = numel(vecAssignPoints)-numel(vecUsePoints);
				if intDeleteUpTo>0
					vecAssignPoints(1:intDeleteUpTo) = []; %remove out-of-bounds entries
				end
			end
			
			%assign data to matrix
			vecUsePoints = intStart:intStop;
			if numel(vecAssignPoints) > numel(vecUsePoints)
				vecAssignPoints = vecAssignPoints(1:numel(vecUsePoints));
			elseif numel(vecAssignPoints) < numel(vecUsePoints)
				vecUsePoints = vecUsePoints(1:numel(vecAssignPoints));
			end
			matTE(intEvent,vecAssignPoints) = vecVals(vecUsePoints);
		else
			vecTheseEdgesT = vecEvents(intEvent) + vecWindow;
			[vecCounts,vecMeans] = makeBins(vecTime,vecVals,vecTheseEdgesT);
			matTE(intEvent,:) = vecMeans;
		end
	end
end

