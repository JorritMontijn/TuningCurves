%% remove all old data
clc
clearvars -except ses
%close all

%% set which analyses to run
sParams.boolDoStimCorrPlot = true;
sParams.boolDoStimDecoding = true;%true;
sParams.boolDoTraceDecoding = true;%true;
sParams.boolDoTracePresenceDecoding = false;
sParams.boolDoTracePEP = true;%true; %peri-event plot
sParams.boolDoStimDecodingPEP = true;%true;
sParams.boolDoOSIPEP = true; %
sParams.boolDoRectify = false;

%% set parameters
sParams.boolMultiPlot = false;

%% set which recordings to process
intProcess=-1;
if intProcess == -1
	%guido's recording
	cellSes{1} = 'Guido';
	cellRec{1} = '20140520xyt04_ses';
elseif intProcess == 0
	cellSes{1} = '20120720';
	cellRec{1} = 'xyt01';
elseif intProcess == 1
	%behavior
	cellSes{1} = '20130627';
	cellSes{2} = '20130627';
	cellRec{1} = 'xyt01';
	cellRec{2} = 'xyt02';
elseif intProcess == 2
	%dir tuning
	%np: 60, 39, 50, 70, 80, 66
	%s: 98, 75, 90, 96, 99, 89
	%s-np: 100, 98, 100, 99, 100, 100
	cellSes{1}= '20120718';
	cellRec{1}= 'xyt01';
	
	cellSes{2}= '20120718';
	cellRec{2}= 'xyt03';
	
	cellSes{3}= '20120720';
	cellRec{3}= 'xyt01';
	
	cellSes{4}= '20120720';
	cellRec{4}= 'xyt03';
	
	cellSes{5}= '20121207';
	cellRec{5}= 'xyt01';
	
	cellSes{6}= '20121207';
	cellRec{6}= 'xyt02';
	%{
	cellSes{7}= '20120720';
	cellRec{7}= 'xyt01';
	
	%}
elseif intProcess == 3
	%ori tuning anesth
	%np: 24, 29
	%s: 80, 63
	%s-np: 99, 95
	%{
	cellSes{1}= '20130612';
	cellRec{1}= 'xyt02';
	
	cellSes{2}= '20130625';
	cellRec{2}= 'xyt02';
	
	cellSes{3} = '20131022';
	cellRec{3} = 'xyt01';
	
	cellSes{4} = '20131022';
	cellRec{4} = 'xyt03';
	%}
	cellSes{1} = '20140129';
	cellRec{1} = 'xyt01';
	
	cellSes{2} = '20140131';
	cellRec{2} = 'xyt01';
elseif intProcess == 4
	%plaids
	cellSes{1}= '20130612';
	cellRec{1}= 'xyt01';
	
	cellSes{2}= '20130625';
	cellRec{2}= 'xyt01';
elseif intProcess == 5
	%old dir/plaids
	cellSes{1}= '20130313';
	cellRec{1}= 'xyt02';
	
	cellSes{2}= '20130313';
	cellRec{2}= 'xyt03';
end

%% get data
sParams.strFigDir = 'D:\Data\Results\neuropilTimeCourses\';
boolAwake = 0;
strAwake = 'anesth';
for boolDoShuffle = [0 1]
	sMultiData = struct;
	
	for intRec=1:numel(cellSes)
		%% load data
		if intProcess >= 0
			clear('ses','sParamsSD','sOutDB','sOutTD','sTC','sBR')
			strSes = cellSes{intRec};
			strRec = cellRec{intRec};

			fprintf('\n   >>> Loading and preparing data file for session %s%s\n\n',strSes,strRec);

			strDir = ['D:' filesep 'Data' filesep 'Processed' filesep 'imagingdata' filesep strSes filesep strRec filesep];
			strFile = sprintf('%s%s_ses.mat',strSes,strRec);

			load([strDir strFile]);
		end
		
		%get shuffling vector
		if boolDoShuffle
			intStims = length(ses.structStim.FrameOn);
			vecShuffle = randperm(intStims);
		else
			vecShuffle = [];
		end
		
		%% perform analysis
		for intDataType=[5 4 3]
			%close figures
			close all;
			drawnow;
			
			%get type string and perform data transformation
			if intDataType == 0
				strDataType='snp';
			elseif intDataType == 3
				strDataType='soma';
				ses = doRecalcdFoF(ses,3); %3=s;4=np;5=s-np
			elseif intDataType == 4
				strDataType='np';
				ses = doRecalcdFoF(ses,4); %3=s;4=np;5=s-np
			elseif intDataType == 5
				strDataType='snp';
				ses = doRecalcdFoF(ses,5); %3=s;4=np;5=s-np
				[ses,indKeepList] = doRecalcdFoF(ses,6);
			end
			
			%set parameters
			sParams.boolUseBaseline = true;
			sParams.intDataType = intDataType;
			sParams.vecShuffle = vecShuffle;
			sParams.vecLikelihoodTrials = 1:round((length(ses.structStim.FrameOn)/2));
			
			%perform analyses
			sProcData = doStimAnalysis(ses,sParams);
			sProcData.indKeepList = indKeepList;
			strSes = ['ses' ses.session];
			strRec = ses.recording;
			if isnumeric(strRec),strRec = sprintf('xyt%02.0f',strRec);end
			sMultiData.([strSes strRec]).(strAwake).(strDataType).sProcData = sProcData;
		end
	end
	if boolDoShuffle
		strShuffled = '_shuffled';
	else
		strShuffled = '_nonshuffled';
	end
	strFilename = sprintf('%sprocessedMultiData%s_trialsSplit%s',strrep(strrep(datestr(clock),':',''),' ','_'),strShuffled,strAwake);
	save(strFilename,'sMultiData','-v7.3')
end