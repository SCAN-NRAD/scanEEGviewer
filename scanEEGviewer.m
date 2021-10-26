
%%%
% scanEEGviewer()							% opens file selection gui
%
% scanEEGviewer(EEGfile)					% absolute or relative path to a .mat file
%												containing at least one variable 'EEG'
%												being a n-by-k matrix where
%												columns correspond to channels and
%												rows correspond to timepoints												
% scanEEGviewer(EEGfile, eventFile)			% absolute or relative paths to two .mat files
%												containing EEG and associated event data
% scanEEGviewer(EEGdata, fs)				% 'EEGdata' is a n-by-k matrix and
%												'fs' is the data's sampling frequency
% scanEEGviewer(EEGdata, fs, labels)		% 'labels' is a cell array of length k
%												containing the EEG channels' names
% scanEEGviewer(EEGdata, fs, labels, eventFile)
%
%
% (c) Michael Müller,	May 2020
%	ORCiD: 0000-0002-6915-4820
%	last modification: 26.10.2021
%%%

function scanEEGviewer(varargin)
	%%% default data path for relative path specifications
	dataPath = './Data/';
	eventPath = 'annot*';
	%%%
	
	%%% get data location
	eventFileFlag = false;
	if(nargin < 1)
		% open file selection gui
		[fileName, dataPath] = uigetfile('*.mat', 'Select an EEG .mat file', dataPath);
		if(~fileName)
			return
		end
	elseif(ischar(varargin{1}))
		fileName = varargin{1};
		% separate string in case of absolute path
		if(strcmp(fileName(1), '/') || strcmp(fileName(1), '.') || strcmp(fileName(1), '..') ...	% linux
				|| strcmp(fileName(2:3), ':\') || strcmp(fileName(2:3), ':/'))						% windows
			endPathIdx = regexp(fileName, {'\', '/'});
			endPathIdx = sort([endPathIdx{:}]);
			dataPath = fileName(1:endPathIdx(end));
			fileName = fileName(endPathIdx(end)+1:end);
		end
	end
	
	%%% load data
	fWB = waitbar(0, 'Loading data...');
	if(exist('fileName', 'var'))
		% if file paths are provided
		try
			data = load([dataPath fileName]);
			rawEEGdata = data.EEG;
		catch
			close(fWB);
			error('Data load error!');
		end
		% read sampling frequeuncy
		try
			fs = data.fs;
		catch
			fs = NaN;
		end
		% read channel labels or use numbering
		try
			chLabel = data.ELabel;
		catch
			% enumerate channels if no labels are provided
			chLabel = cellstr(num2str((1:size(rawEEGdata, 2))'));
		end
		% try to read resected channels from data file
		try
			RBTchLab = data.RBTchLab;
			[~, RBTchIdx] = intersect(chLabel, RBTchLab);
			if(numel(RBTchIdx) < numel(RBTchLab))
				warndlg('Not all RBT channels could be identified.', 'Data read error');
			end
		catch
			RBTchIdx = NaN;
		end
		
		if(nargin > 1)
			% set path to event file if specified
			eventFilename = varargin{2};
			eventFileFlag = true; % explicitely set event file
		else
			% use default path to event file to check
			eventFilename = [fileName(1:end-4) '-events*'];
		end
	else
		% if data is provided
		fileName = [datestr(date(), 'yymmdd') '.mat'];	% set dummy file name
		rawEEGdata = varargin{1};
		RBTchIdx = NaN;
		
		if(nargin < 2)
			fs = NaN;
		else
			fs = varargin{2};
		end
		if(nargin < 3 || isempty(varargin{3}))
			% enumerate channels if no labels are provided
			chLabel = cellstr(num2str((1:size(rawEEGdata, 2))'));
		else
			chLabel = varargin{3};
		end
		if(nargin > 3)
			% set path to event file if specified
			eventFilename = varargin{4};
			eventFileFlag = true;
		else
			eventFilename = [];
		end
	end
	
	% if fs is unknown request input
	while(isnan(fs))
		tmpFs = inputdlg('Please specify the data''s sampling frequency');
		if(isempty(tmpFs))
			close(fWB);
			return
		end
		fs = str2double(tmpFs{1});
	end
	
	[nbSamp, nbChan] = size(rawEEGdata);
	
	% check data consistency and apply simple fixes if necessary
	if(nbChan < numel(chLabel))
		warning(['Number of EEG channels and labels not consistent! ' ...
				'Labels ' num2str(nbChan+1) ' to ' num2str(numel(chLabel)) ...
				' (' [chLabel{nbChan+1:end}] ') ignored.']);
		chLabel = chLabel(1:nbChan);
	elseif(numel(chLabel) < nbChan)
		warning(['Number of EEG channels and labels not consistent! ' ...
				'Labels ' num2str(numel(chLabel)+1) ' to ' num2str(nbChan) ' padded.']);
		chLabel(end+1:nbChan) = cellstr(num2str(((numel(chLabel)+1):nbChan)'));
	end
% 	labHyph = cellfun(@(x) strfind(x, '-'), chLabel, 'uni',false);
% 	labUndl = cellfun(@(x) strfind(x, '_'), chLabel, 'uni',false);
% 	if(any(diff(cellfun(@numel, labHyph))) || any(diff(cellfun(@numel, labUndl))))
% 		warning('Style of EEG labels is not consistent! Using enumerated labels instead.');
% 		chLabel = cellstr(num2str((1:size(rawEEGdata, 2))'));
% 	end
	%%%
	
	% set default values
	offset = 300;				% signal display separation
	winLength = 10;				% window length [s]
	winStartIdx = 0;			% window start value [s]
	ovWinLen = 1;				% win length in spectral overview [m]
	ovWinLines = cell(2, 1);
	
	% set initial time slider values
	sldMin = 0;
	sldMax = nbSamp / fs;
	sldStep = [1 winLength] * fs / nbSamp; % short & long step
	
	% set initial data to display
	filtEEGdata = rawEEGdata;
	loCut = 0;	% filter
	hiCut = 0;	% filter
	dispEEGdata = rawEEGdata;
	dispChan = 1:nbChan;
	dispChLabel = chLabel;
% 	channelLat = zeros(numel(chLabel), 1);
	
	% display colors
	elCols = [0 0 1; 1 0 0; 0 1 0; 1 0 1; 1 1 0; 0 1 1];	% signals
	boxCol = {'r', '', 'm'};								% boxes
	
	% allocate 
	[elUni, elGroup] = deal([]);
	[specCA, specBI] = deal([]);	% storing spectrum data of both referencing types
	biChCombs = [];
	boxCoord = [];
	boxObj = [];
	boxBtn = 0;
	
	% setup event list	
	eventTypes = {'Spike(s)', 'Slow wave(s)', 'Sharp wave(s)', ...
					'Fast oscillations', 'Analysis end', 'General note'};	% default event types
	eventList = {};
	selEvType = true(1, 6);
	
	%%% setup GUI
	fEvent = figure('Visible','off', 'numbertitle','off');
	fMain = figure('Visible','off', 'Units','normalized','Position',[0.01, 0.05, 0.98, 0.85], ...
					'Name',['SCAN EEG viewer - ' fileName], 'numbertitle','off', 'KeyPressFcn',@keyDownListener);
	fMain.WindowButtonMotionFcn = @mouseDown;
	
	axEEG = axes('Position',[0.04 0.07 0.79 0.8], 'ClippingStyle','rectangle');
	axEEG.ButtonDownFcn = @signalMouse;
	
	axSpec = axes('Position',[0.04 0.90 0.79 0.09]);
	axSpec.ButtonDownFcn = @specGoto;
	
	cbLogSpec		= uicontrol('Parent',fMain, 'Style','checkbox', 'String','log.', ...
								'Units','normalized', 'Position',[0.0025 0.93 0.025 0.0375], ...
								'Callback',@switchSpecScale);
	
	cbColGrad		= uicontrol('Parent',fMain, 'Style','checkbox', 'String','col.grad.', ...
								'Units','normalized', 'Position',[0.0025 0.87 0.05 0.015], ...
								'Value',1, 'Callback',@switchDispOpt);
	
	tbFilt			= uicontrol('Parent',fMain, 'Style','togglebutton', 'String','Filter', ...
								'Units','normalized', 'Position',[0.85 0.95 0.04 0.03], ...
								'Callback',@switchFilter);
	edtFiltLC		= uicontrol('Parent',fMain, 'Style','edit', 'String','0.5', ...
								'Units','normalized', 'Position',[0.90 0.95 0.04 0.03], ...
								'Callback',@setFilterCut);
	edtFiltHC		= uicontrol('Parent',fMain, 'Style','edit', 'String','50', ...
								'Units','normalized', 'Position',[0.95 0.95 0.04 0.03], ...
								'Callback',@setFilterCut);
	
	tbRefGroup		= uibuttongroup('Parent',fMain, 'Units','normalized', 'Position',[0.85 0.9 0.14 0.04]);
	tbRefRec		= uicontrol(tbRefGroup, 'Style','togglebutton', 'String','Recording', ...
								'Units','normalized', 'Position',[0.01 0.05 0.32 0.9], ...
								'Callback',@changeReference);
	tbRefCA			= uicontrol(tbRefGroup, 'Style','togglebutton', 'String','Common Av.', ...
								'Units','normalized', 'Position',[0.34 0.05 0.32 0.9], ...
								'Callback',@changeReference);
	tbRefBI			= uicontrol(tbRefGroup, 'Style','togglebutton', 'String','Bipolar', ...
								'Units','normalized', 'Position',[0.67 0.05 0.32 0.9], ...
								'Callback',@changeReference);	
	
	btnSensD		= uicontrol('Parent',fMain, 'Style','pushbutton', 'String','-', ...
								'Units','normalized', 'Position',[0.85 0.84 0.065 0.03], ...
								'Callback',@changeSens);
	btnSensU		= uicontrol('Parent',fMain, 'Style','pushbutton', 'String','+', ...
								'Units','normalized', 'Position',[0.925 0.84 0.065 0.03], ...
								'Callback',@changeSens);
	
    lbxEEGchan		= uicontrol('Parent',fMain, 'Style','listbox', ...
								'String',chLabel, ...
								'Min',1,'Max',nbChan,'Value',1:nbChan, ...
								'Units','normalized', 'Position',[0.85 0.07 0.04 0.75], ...
								'Callback',@selEEGchan);
	
	txtEvents		= uicontrol('Parent',fMain, 'Style','text', 'String',[], ...
								'Units','normalized', 'Position',[0.90 0.78 0.09 0.04], ...
								'FontWeight','bold');
	cbSpike			= uicontrol('Parent',fMain, 'Style','checkbox', 'String',eventTypes{1}, ...
								'Units','normalized', 'Position',[0.90 0.765 0.09 0.015], ...
								'Value',1, 'Callback',@switchEvtSel);
	cbSlowW			= uicontrol('Parent',fMain, 'Style','checkbox', 'String',eventTypes{2}, ...
								'Units','normalized', 'Position',[0.90 0.75 0.09 0.015], ...
								'Value',1, 'Callback',@switchEvtSel);
	cbSharpW		= uicontrol('Parent',fMain, 'Style','checkbox', 'String',eventTypes{3}, ...
								'Units','normalized', 'Position',[0.90 0.735 0.09 0.015], ...
								'Value',1, 'Callback',@switchEvtSel);
	cbFastOsc		= uicontrol('Parent',fMain, 'Style','checkbox', 'String',eventTypes{4}, ...
								'Units','normalized', 'Position',[0.90 0.720 0.09 0.015], ...
								'Value',1, 'Callback',@switchEvtSel);
	lbxEvents		= uicontrol('Parent',fMain, 'Style','listbox', ...
								'String',eventList, ...
								'Min',0,'Max',1,'Value',0, ...
								'Units','normalized', 'Position',[0.90 0.07 0.09 0.64], ...
								'Callback',@gotoEvent);
	btnEvDel		= uicontrol('Parent',fMain, 'Style','pushbutton', ...
								'String','Delete event', ...
								'Units','normalized', 'Position',[0.94 0.03 0.05 0.03], ...
								'Callback',@delEvent);
	
	edtwinL			= uicontrol('Parent',fMain, 'Style','edit', 'String',winLength, ...
								'Units','normalized', 'Position',[0.0075 0.01 0.025 0.03], ...
								'Callback',@changeWinLen);
	sldTime			= uicontrol('Parent',fMain, 'Style','slider', ...
								'Min',sldMin,'Max',sldMax,'Value',sldMin,'SliderStep',sldStep, ...
								'Units','normalized','Position',[0.04 0.01 0.79 0.03], ...
								'Callback',@selTime);
	
	cbSOZ			= uicontrol('Parent',fMain, 'Style','checkbox', 'String','SOZ', ...
								'Units','normalized', 'Position',[0.85 0.025 0.09 0.015], ...
								'Value',0, 'Enable','off', 'Callback',@switchDispOpt);
	cbRBT			= uicontrol('Parent',fMain, 'Style','checkbox', 'String','RBT', ...
								'Units','normalized', 'Position',[0.85 0.01 0.09 0.015], ...
								'Value',0, 'Callback',@switchDispOpt);
	if(isnan(RBTchIdx))
		set(cbRBT, 'Enable','off');
	end
	
	% set events/SOZ if present
	raterID = [];	% ID of person making annotations
	SOZchIdx = [];
	try
		if(~eventFileFlag)
			eventDir = dir([dataPath eventPath]);
			eventPath = [];

			if(numel(eventDir) > 0 && ~isempty(eventFilename))

				% check all annotation directories for valid event files
				eventFiles = cellfun(@(x) dir([dataPath x '/' eventFilename]), {eventDir.name}, 'uni',false);

				allEvtFiles = cell(0);
				for dirIdx = 1:numel(eventDir)
					for fIdx = 1:numel(eventFiles{dirIdx})
						allEvtFiles{end+1} = [eventDir(dirIdx).name '/' eventFiles{dirIdx}(fIdx).name];
					end
				end

				[selIdx, selFlag] = listdlg('Name','Event selection', 'ListString',allEvtFiles, ...
					'PromptString','Select event file', 'SelectionMode','single', 'ListSize',[350 150]);

				if(selFlag)
					eventFileFlag = true;

					% load data
					eventFilename = strsplit(allEvtFiles{selIdx}, '/');
					eventPath = [eventFilename{1} '/'];
					eventFilename = eventFilename{2};
					eventData = load([dataPath eventPath eventFilename]);
				end
			end
		else
			eventData = load(eventFilename);
		end
		
		if(~strcmp(eventData.fileName, fileName))
			warndlg('Event data filename does not match data filename.', 'Data missmatch');
		end

		raterID = eventData.raterID;
		eventList = eventData.eventList;
		eventTypes = eventData.eventTypes;

		% update GUI elements
		set(fMain, 'Name',[get(fMain, 'Name') ' - annotations by ' raterID]);
		setEventfileNameTxt(33);
		updateEventListbox(numel(eventList));

		% try to read SOZ channels from 'End Analysis' event comment
		try
			endAnaIdx = cellfun(@(x) x{1}(5), eventList);
			SOZchName = strtrim(strsplit(eventList{endAnaIdx}{5}, ','));
			[~, SOZchIdx] = intersect(chLabel, SOZchName);
			if(~isempty(SOZchIdx))
				set(cbSOZ, 'Enable','on');
			end
		catch
			% simply here to please the machine
		end
	catch
		if(eventFileFlag)
			msgbox('Event file could not be loaded.', 'Error','error');
		end
	end
	
	waitbar(0.25, fWB, 'Filtering data...');
	filterEEGsignals();
	
	% pre-calculate spectra for recording & bipolar reference
	calcSpectra();

	close(fWB);
	set(fMain, 'Visible','on');
	
	
	%%% callback functions
	function switchSpecScale(~, ~)
		% switch spectrum color-scale lin/log
		plotSpectrum();
	end
	function switchDispOpt(~, ~)
		% switch electrode-wise signal colors gradient/static
		plotEEGsignals();
	end
	function switchFilter(~, ~)
		% switch signal bandpass filter on/off
		if(get(tbFilt, 'Value'))
			filterEEGsignals();
		end
		changeEEGsigDisp();
	end
	function setFilterCut(~, ~)
		% set bandpass filter values and filter on
		set(tbFilt, 'Value',1);
		switchFilter();
	end
	function selEEGchan(~, ~)
		% set selected channels to display and adapt event list display
		changeEEGsigDisp();
		updateEventListbox(get(lbxEvents, 'Value'));
	end
	function changeReference(~, ~)
		% set signal reference
		changeEEGsigDisp();
	end
	function selTime(source, ~)
		% set time window of signal display
		winStartIdx = get(source, 'Value');
		plotSpectrumPosition();
		plotEEGsignals();
	end
	function switchEvtSel(~, ~)
		% set event types to display
		selEvType = logical([get(cbSpike, 'Value') ...
								get(cbSlowW, 'Value') ...
								get(cbSharpW, 'Value') ...
								get(cbFastOsc, 'Value') ...
								true true]);
		plotEEGsignals();
		updateEventListbox(get(lbxEvents, 'Value'));
	end
	function changeWinLen(source, ~)
		% change signal window length
		newVal = str2double(get(source, 'String'));
		if(isnan(newVal))
			% reset previous value
			set(edtwinL, 'String',winLength);
			return
		end
		
		winLength = newVal;
		set(sldTime, 'SliderStep',[winLength/10 winLength] * fs / nbSamp);
		plotSpectrumPosition();
		plotEEGsignals();
	end
	function changeSens(source, ~)
		% change signal amplitude scale
		switch(source.String)
			case '+'
				offset = offset / 1.5;
			case '-'
				offset = offset * 1.5;
		end
		plotEEGsignals();
	end
	function signalMouse(~, event)
		if(isvalid(fEvent))
			if(isempty(boxObj))
				boxObj = rectangle();	% add dummy object
			end
			close(fEvent);
		end
		
		clickDown = get(axEEG, 'CurrentPoint');
		clickDown = clickDown(1, 1:2);
		
		if(clickDown(1) > sldMax)
			return
		end
		
		if(~isempty(boxCoord))
			boxBtn = 0;
			if(event.Button == 1)		% left click
				boxCoord = [];
			elseif(event.Button == 3)	% right click
				% get in-boxed time
				timePos = sort([boxCoord(1) clickDown(1)]);
				timePos = round(timePos * 100) / 100;
				
				% get in-boxed channels
				selCh = sort([boxCoord(2) clickDown(2)]);
				selCh = max(ceil(selCh(1)/offset), 1):min(floor(selCh(2)/offset), size(dispEEGdata, 2));
				if(isempty(selCh))
					selCh = min(max(round(mean([boxCoord(2) clickDown(2)])/offset), 1), size(dispEEGdata, 2));
				end
				
				if(tbRefBI.Value)
					selCh = biChCombs(selCh, :);
					selCh = unique(selCh);
				end
				
				boxCoord = [];
				eventWin(-1, timePos(1), timePos(2), selCh);
			end
		elseif(~isempty(boxObj))
			delete(boxObj);
			boxObj = [];
		else
			boxBtn = event.Button;
			boxCoord = clickDown;
		end
	end
	function mouseDown(~, ~)
		if(~isempty(boxCoord))
			delete(boxObj);
			
			currPoint = get(axEEG, 'CurrentPoint');
			currPoint = currPoint(1, 1:2);		
			
			boxDim = abs(currPoint - boxCoord);
			rCoord = min([boxCoord; currPoint]);
			
			switch(boxBtn)
				case 1
					% visual aid box
					yLab = [num2str(round(boxDim(2))) ' \muV'];
				case 3
					% capture event box
					selCh = sort([boxCoord(2) currPoint(2)]);
					firstCh = ceil(selCh(1) / offset);
					lastCh = floor(selCh(2) / offset);
					
					if(lastCh < firstCh)
						chIdx = min(max(round(mean(selCh) / offset), 1), size(dispEEGdata, 2));
					else
						chIdx = unique(min(max([firstCh lastCh], 1), size(dispEEGdata, 2)));
					end
					
					if(tbRefBI.Value)
						chIdx = biChCombs(chIdx, :);
						chIdx = unique(chIdx);
						chIdx = chIdx([1 end]);
					end
					
					yLab = strjoin(chLabel(dispChan(chIdx)), ' - ');
			end
			
			r = rectangle('Position',[rCoord boxDim], ...
							'EdgeColor',boxCol{boxBtn}, 'FaceColor','none', ...
							'Parent',axEEG, 'HitTest','off');
			
			if(rCoord(1) > winStartIdx && ...
					rCoord(1) + boxDim(1) > winStartIdx + winLength * 0.95)
				xVal = rCoord(1) - winLength / 20;
			else
				xVal = rCoord(1) + boxDim(1) + winLength/100;
			end
			ySel = text(xVal, rCoord(2) + boxDim(2)/2, yLab, ...
						'Color',boxCol{boxBtn}, 'BackgroundColor',[1 1 1 0.67], ...
						'FontSize',12, 'FontWeight','bold', 'Parent',axEEG, 'HitTest','off');
			
			winH = get(axEEG, 'YLim');
			xMid = (max(rCoord(1), winStartIdx) + min(rCoord(1)+boxDim(1), winStartIdx+winLength)) / 2;
			tDur = text(xMid - winLength/100, rCoord(2) - winH(2)/50, ...
						[num2str(round(boxDim(1) * 1000)) ' ms'], ...
						'Color',boxCol{boxBtn}, 'BackgroundColor',[1 1 1 0.67], ...
						'FontSize',12, 'FontWeight','bold', 'Parent',axEEG, 'HitTest','off');
			
			boxObj = [r, ySel, tDur];
		end
	end


	%%% EEG signal functions
	function filterEEGsignals()
		newLoCut = str2double(get(edtFiltLC, 'String'));
		newHiCut = str2double(get(edtFiltHC, 'String'));
		
		if(newLoCut ~= loCut || newHiCut ~= hiCut)
			loCut = newLoCut;
			hiCut = newHiCut;
			[B, A] = butter(4, [loCut hiCut] * 2 / fs);
			filtEEGdata = filtfilt(B, A, rawEEGdata);
		end
	end
	function changeEEGsigDisp()
		if(isvalid(fEvent))
			close(fEvent);
		end
		
		if(get(tbFilt, 'Value'))
			tmpEEGdata = filtEEGdata;
		else
			tmpEEGdata = rawEEGdata;
		end
		
		dispChan = get(lbxEEGchan, 'Value');
		dispChLabel = chLabel(dispChan);

		% use label letters to discriminate electrodes
		labelVal = cellfun(@(x) double(lower(x)), dispChLabel, 'uni',false);
		labelChar = cellfun(@(x) x(x > 96 & x < 123), labelVal, 'uni',false);

		labelChar = cellfun(@sum, labelChar);
		[elUni, ~, elGroup] = unique(labelChar, 'stable');
		
		% set signal reference
		if(get(tbRefRec, 'Value'))
			% use record reference
			dispEEGdata = tmpEEGdata(:, dispChan);			
		elseif(get(tbRefCA, 'Value'))
			% use common average of selected channels only
			dispEEGdata = bsxfun(@minus, tmpEEGdata(:, dispChan), mean(tmpEEGdata(:, dispChan), 2));
		elseif(get(tbRefBI, 'Value'))
			% use bipolar referencing
			labelCharSwitch = diff(labelChar);
			
			biChIdx = find(labelCharSwitch == 0);
			biChCombs = [biChIdx+1 biChIdx];
			
			dispEEGdata = tmpEEGdata(:, dispChan(biChIdx+1)) - tmpEEGdata(:, dispChan(biChIdx));
			dispChLabel = cellfun(@(x, y) [x '-' y], dispChLabel(biChIdx+1), ...
									dispChLabel(biChIdx), 'uni',false);
% 			channelLat = channelLat(biChIdx);
			[elUni, ~, elGroup] = unique(labelChar(biChIdx), 'stable');
		end
		
		plotSpectrum();
		plotEEGsignals();
	end
	function plotEEGsignals()
		winEndIdx = winStartIdx + winLength;
		selDataIdx = max(winStartIdx*fs, 1):min(winEndIdx*fs+1, nbSamp);
		
		if(isempty(dispChLabel))
			return
		end
		
		nbDispCh = numel(dispChLabel);
		
		tmpEEGdata = dispEEGdata(round(selDataIdx), :);
		tmpEEGdata = bsxfun(@plus, (1:nbDispCh)*offset, tmpEEGdata);
		
		cla(axEEG);
		hold(axEEG, 'on');
		
		% vertical lines every second
		lineXpos = round(winStartIdx:winEndIdx);
		for i = 1:numel(lineXpos)
			line([lineXpos(i) lineXpos(i)], [-offset/4 (nbDispCh+1)*offset], ...
					'Color',[0.75 0.75 0.75], 'LineStyle','--', ...
					'Parent',axEEG, 'HitTest','off');
		end
		
		nbEl = numel(elUni);
		colIdx = mod((1:nbEl)-1, size(elCols, 1))+1;
		
		if(cbColGrad.Value)
			% plot channels with color gradient within electrodes
			for elIdx = 1:nbEl
				elCh = find(elGroup == elIdx);
				elNbCh = numel(elCh);
				chColShades = fliplr(linspace(max(0.25, 0.75-0.1*elNbCh), 0.75, elNbCh));

				for chIdx = 1:elNbCh
					plot(axEEG, selDataIdx/fs, tmpEEGdata(:, elCh(chIdx)), ...
							'Color',elCols(colIdx(elIdx), :)*chColShades(chIdx), 'HitTest','off');
				end
			end
		else
			% plot entire electrodes in one color (faster)
			for elIdx = 1:nbEl
				plot(axEEG, selDataIdx/fs, tmpEEGdata(:, elGroup == elIdx), ...
						'Color',elCols(colIdx(elIdx), :)*0.5, 'HitTest','off');
			end
		end
		
		% mark channel labels if RBT/SOZ
		markLabel = dispChLabel;
		if(~get(tbRefBI, 'Value'))
			if(get(cbRBT, 'Value'))
				[~, rbtIdx] = intersect(dispChan, RBTchIdx);
				markLabel(rbtIdx) = cellfun(@(x) ['\bf{' x '}_{RBT}'], markLabel(rbtIdx), 'uni',false);
			end
			if(get(cbSOZ, 'Value'))
				[~, sozIdx] = intersect(dispChan, SOZchIdx);
				markLabel(sozIdx) = cellfun(@(x) [x '^{SOZ}'], markLabel(sozIdx), 'uni',false);
			end
		else
			set(cbSOZ, 'Value',0);
			set(cbRBT, 'Value',0);
		end
		
		% set axis
		set(axEEG, 'YDir','reverse', 'YTick',(1:nbDispCh)*offset, 'YTickLabel',markLabel, 'FontSize',8, 'box','on');
		xlim(axEEG, [winStartIdx winEndIdx]);
		ylim(axEEG, [-offset/4 (nbDispCh+1) * offset]);
		
		yLimits = get(axEEG, 'YLim');
		
		% plot time meter
		xPos = winStartIdx + winLength * 0.025;
		yPos = yLimits(2) * [0.93 0.97];
		line([xPos xPos+0.07], [yPos(1) yPos(1)], 'Color','k', 'LineWidth',3, 'Parent',axEEG, 'HitTest','off');
		line([xPos xPos+0.2], [yPos(2) yPos(2)], 'Color','k', 'LineWidth',3, 'Parent',axEEG, 'HitTest','off');
		text(xPos + 0.07 + winLength/200, yPos(1), '70 ms', 'FontSize',11, 'FontWeight','bold', ...
			'BackgroundColor',[1 1 1 0.67], 'Parent',axEEG, 'HitTest','off');
		text(xPos + 0.2 + winLength/200, yPos(2), '200 ms', 'FontSize',11, 'FontWeight','bold', ...
			'BackgroundColor',[1 1 1 0.67], 'Parent',axEEG, 'HitTest','off');
		
		% plot amplitude meter
		xPos = winStartIdx + winLength * [0.935 0.95];
		yPos = yLimits(2)*0.95 + [-300 0];
		line(xPos, [yPos(1) yPos(1)], 'Color','k', 'LineWidth',2, 'Parent',axEEG, 'HitTest','off');
		line(xPos, [yPos(2) yPos(2)], 'Color','k', 'LineWidth',2, 'Parent',axEEG, 'HitTest','off');
		line(mean(xPos) * [1 1], [yPos(1) yPos(2)], 'Color','k', 'LineWidth',2, 'Parent',axEEG, 'HitTest','off');
		text(xPos(1) - winLength/37.5, yPos(2), '0 \muV', 'FontSize',11, 'FontWeight','bold', ...
			'BackgroundColor',[1 1 1 0.67], 'Parent',axEEG, 'HitTest','off');
		text(xPos(2) + winLength/200, yPos(1), '300 \muV', 'FontSize',11, 'FontWeight','bold', ...
			'BackgroundColor',[1 1 1 0.67], 'Parent',axEEG, 'HitTest','off');
		
		if(isempty(boxCoord))
			boxObj = [];
		else
			mouseDown([], []);
		end
		plotEvents();
	end


	%%% spectrum functions
	function calcSpectra()
		specFreq = 0:0.5:30;
		nbWin = floor(nbSamp / (fs * ovWinLen));
		
		% set bipolar reference and calulate spectrum
		waitbar(0.50, fWB, 'Calculate spectra...');
		
		set(tbRefBI, 'Value',1);
		changeEEGsigDisp();
		
		tmpSpec = deal(zeros(numel(specFreq), nbWin, size(dispEEGdata, 2)));
		for chIdx = 1:size(dispEEGdata, 2)
			s = spectrogram(dispEEGdata(:, chIdx), fs * ovWinLen, 0, specFreq, fs);
			tmpSpec(:, :, chIdx) = abs(s);
		end
		specBI = mean(tmpSpec, 3);
		
		% calculate recorde reference spectrum and display it
		waitbar(0.75, fWB, 'Calculate spectra...');
		
		tmpSpec = deal(zeros(numel(specFreq), nbWin, size(rawEEGdata, 2)));
		for chIdx = 1:size(rawEEGdata, 2)
			s = spectrogram(rawEEGdata(:, chIdx), fs * ovWinLen, 0, specFreq, fs);
			tmpSpec(:, :, chIdx) = abs(s);
		end
		specCA = mean(tmpSpec, 3);
		
		set(tbRefRec, 'Value',1);
		changeEEGsigDisp();
	end
	function plotSpectrum()
		% plot spectrum
		if(get(tbRefBI, 'Value'))
			tmpSpec = specBI;
		else
			tmpSpec = specCA;
		end
		
		if(cbLogSpec.Value)
			tmpSpec = log(tmpSpec);
		end
		
		nbWin = floor(nbSamp / (fs * ovWinLen));
		
		imagesc((0:(nbWin-1) * ovWinLen)/60, 0:0.5:30, tmpSpec, 'Parent',axSpec, 'HitTest','off');
		colormap(jet);
		axis(axSpec, 'tight');
		set(axSpec, 'YDir','normal', 'ButtonDownFcn',@specGoto, 'box','on');
		
		plotSpectrumPosition();
	end
	function plotSpectrumPosition()
		% plot vertical lines on spectrum indicating current position
		delete(ovWinLines{1});
		delete(ovWinLines{2});
		
		xLineStart = (winStartIdx - ovWinLen/2) / 60;
		xLineEnd = (winStartIdx + winLength - ovWinLen/2) / 60;
		
		yLimVal = get(axSpec, 'YLim');
		l1 = line([xLineStart xLineStart], yLimVal, 'Color','r', 'LineWidth',2, 'Parent',axSpec, 'HitTest','off');
		l2 = line([xLineEnd xLineEnd], yLimVal, 'Color','r', 'LineWidth',2, 'Parent',axSpec, 'HitTest','off');
		
		ovWinLines = {l1, l2};
	end
	function specGoto(~, event)
		% change current position to center clicked point on spectrum
		timePos = event.IntersectionPoint(1) * 60 - winLength/2;
		timePos = max(timePos, 0);
		winStartIdx = min(timePos, sldMax);
		set(sldTime, 'Value',winStartIdx);
		plotSpectrumPosition();
		plotEEGsignals();
	end


	%%% event functions
	function eventWin(evIdx, evStart, evEnd, selCh)
		% request annotater ID if not yet defined
		if(isempty(raterID))
			tmpID = inputdlg('Please enter an ID:');
			if(isempty(tmpID) || isempty(tmpID{1}))
				return
			end
			raterID = tmpID{1};
			set(fMain, 'Name',[get(fMain, 'Name') ' - annotations by ' raterID]);
			
			eventPath = ['annot_' raterID '/'];
			if(~exist([dataPath eventPath], 'dir'))
				mkdir([dataPath eventPath]);
			end
			eventFilename = [fileName(1:end-4) '-events-' raterID];
			
			multFcount = 1;
			tmpName = eventFilename;
			while(~isempty(dir([dataPath eventPath tmpName '.mat'])))
				multFcount = multFcount + 1;
				tmpName = [eventFilename '(' num2str(multFcount) ')'];
			end
			eventFilename = [tmpName '.mat'];
			
			setEventfileNameTxt(33);
		end
		
		if(evIdx > 0)
			winTitle = 'Modify event';
			selEv = eventList{evIdx};
			evType = selEv{1};
			evStart = selEv{2};
			evEnd = selEv{3};
			evDur = evEnd - evStart;
			evChan = selEv{4};
			evComm = selEv{5};
		else
			winTitle = 'Add event';
			evType = false(numel(eventTypes), 1);
			evDur = evEnd - evStart;
			evChan = dispChan(selCh);
			evComm = [];
		end
		
		if(numel(evType) < numel(eventTypes))
			evType = [evType; false(numel(eventTypes) - numel(evType))];
		end
		
		% display event gui near clicked point
		mousePos = get(0, 'PointerLocation');
		mousePos(2) = max(mousePos(2), 500);
		
		fEvent = figure('MenuBar','None', 'NumberTitle','off', 'Name',winTitle, ...
						'Position',[mousePos(1)+50 mousePos(2)-450 300 450]);
		
		cbEvt1		= uicontrol('Parent',fEvent, 'Style','checkbox', 'String',eventTypes{1}, ...
								'Units','normalized', 'Position',[0.03 0.94 0.45 0.03], 'Value',evType(1));
		cbEvt2		= uicontrol('Parent',fEvent, 'Style','checkbox', 'String',eventTypes{2}, ...
								'Units','normalized', 'Position',[0.03 0.90 0.45 0.03], 'Value',evType(2));
		cbEvt3		= uicontrol('Parent',fEvent, 'Style','checkbox', 'String',eventTypes{3}, ...
								'Units','normalized', 'Position',[0.03 0.86 0.45 0.03], 'Value',evType(3));
		cbEvt4		= uicontrol('Parent',fEvent, 'Style','checkbox', 'String',eventTypes{4}, ...
								'Units','normalized', 'Position',[0.03 0.82 0.45 0.03], 'Value',evType(4));
		cbEvt5		= uicontrol('Parent',fEvent, 'Style','checkbox', 'String',eventTypes{5}, ...
								'Units','normalized', 'Position',[0.03 0.78 0.45 0.03], 'Value',evType(5));
		cbEvt6		= uicontrol('Parent',fEvent, 'Style','checkbox', 'String',eventTypes{6}, ...
								'Units','normalized', 'Position',[0.03 0.74 0.45 0.03], 'Value',evType(6));
		
		txtStart	= uicontrol('Parent',fEvent, 'Style','text', 'String','Time point [s]', ...
								'Units','normalized', 'Position',[0.03 0.65 0.33 0.05]);
		edtStart	= uicontrol('Parent',fEvent, 'Style','edit', 'String',evStart, ...
								'Units','normalized', 'Position',[0.38 0.65 0.59 0.05], ...
								'Callback',@editEventStart);
		txtDur		= uicontrol('Parent',fEvent, 'Style','text', 'String','Duration [s]', ...
								'Units','normalized', 'Position',[0.03 0.59 0.33 0.05]);
		edtDur		= uicontrol('Parent',fEvent, 'Style','edit', 'String',evDur, ...
								'Units','normalized', 'Position',[0.38 0.59 0.59 0.05], ...
								'Callback',@editEventDur);
		txtEnd		= uicontrol('Parent',fEvent, 'Style','text', 'String','End point [s]', ...
								'Units','normalized', 'Position',[0.03 0.53 0.33 0.05]);
		edtEnd		= uicontrol('Parent',fEvent, 'Style','edit', 'String',evEnd, ...
								'Units','normalized', 'Position',[0.38 0.53 0.59 0.05], ...
								'Callback',@editEventEnd);
		
		lbxEvtCh	= uicontrol('Parent',fEvent, 'Style','listbox', ...
								'String',chLabel, ...
								'Min',1, 'Max',nbChan, 'Value',evChan, ...
								'Units','normalized', 'Position',[0.03 0.21 0.94 0.28]);
		set(lbxEvtCh, 'ListBoxTop',2); % fix so that next line works in all cases
		set(lbxEvtCh, 'ListBoxTop',max(evChan(1)-1, 1));
		
		edtComm		= uicontrol('Parent',fEvent, 'Style','edit', 'Min',0, 'Max',2, ...
								'HorizontalAlignment','left', 'String',evComm, ...
								'Units','normalized', 'Position',[0.03 0.03 0.94 0.17]);
		
		if(numel(raterID) > 9)
			dispAID = ['ID: ' raterID(1:8) '...'];
		else
			dispAID = ['ID: ' raterID];
		end
		txtAID		= uicontrol('Parent',fEvent, 'Style','text', 'String',dispAID, ...
								'Units','normalized', 'Position',[0.52 0.90 0.45 0.07], ...
								'FontWeight','bold');
		btnStrEv	= uicontrol('Parent',fEvent, 'Style','pushbutton', 'String','Store event', ...
								'Units','normalized', 'Position',[0.52 0.82 0.45 0.07], ...
								'Callback',@storeEvent);
		btnCanEv	= uicontrol('Parent',fEvent, 'Style','pushbutton', 'String','Cancel', ...
								'Units','normalized', 'Position',[0.52 0.74 0.45 0.07], ...
								'Callback',@cancelEvent);
		
		function storeEvent(~, ~)
			% compose event data
			eventType = logical([	cbEvt1.Value; ...
									cbEvt2.Value; ...
									cbEvt3.Value; ...
									cbEvt4.Value; ...
									cbEvt5.Value; ...
									cbEvt6.Value]);
			if(~any(eventType))
				msgbox('Please select an event type!');
				return
			end
			
			eventStart = str2double(edtStart.String);
			eventEnd = str2double(edtEnd.String);
			if(isnan(eventStart) || isnan(eventEnd))
				msgbox('Please select valid time points!');
				return
			end
			
			eventChan = lbxEvtCh.Value;
			if(isempty(eventChan))
				msgbox('Please select a channel!');
				return
			end
			
			eventComm = edtComm.String;
			
			eventObj = {eventType, eventStart, eventEnd, eventChan, eventComm};
			
			if(evIdx > 0)
				eventList{evIdx} = eventObj;
			else
				% add new event
				eventList{end+1} = eventObj;
				evIdx = numel(eventList);
			end
			
			% store event in time-sorted order
			sortIdx = saveEventFile();
			
			% close window and refresh display
			cancelEvent();
			updateEventListbox(find(sortIdx == evIdx));
			plotEEGsignals();
		end
		function editEventStart(source, ~)
			% adapt duration (& end) based on start time provided
			startTime = round(str2double(source.String) * 100) / 100;
			endTime = str2double(get(edtEnd, 'String'));
			if(startTime > endTime)
				endTime = startTime;
				set(edtEnd, 'String',endTime);
			end
			set(edtStart, 'String',startTime);
			set(edtDur, 'String',endTime - startTime);
		end
		function editEventDur(source, ~)
			duration = round(str2double(source.String) * 100) / 100;
			startTime = str2double(get(edtStart, 'String'));
			endTime = str2double(get(edtEnd, 'String'));
			if(duration < 0)
				% adapt start time if provided duration is negative
				set(edtStart, 'String',endTime + duration);
				set(edtDur, 'String',-duration);
			else
				% adapt end time if provided duration is positive
				set(edtEnd, 'String',startTime + duration);
				set(edtDur, 'String',duration);
			end
		end
		function editEventEnd(source, ~)
			% adapt duration (& start) based on end time provided
			endTime = round(str2double(source.String) * 100) / 100;
			startTime = str2double(get(edtStart, 'String'));
			if(endTime < startTime)
				endTime = startTime;
				set(edtEnd, 'String',endTime);
			end
			set(edtEnd, 'String',endTime);
			set(edtDur, 'String',endTime - startTime);	
		end
	end
	function sortIdx = saveEventFile()
		% sort events by time
		allEvStarts = cellfun(@(x) x{2}, eventList);
		[~, sortIdx] = sort(allEvStarts);
		eventList = eventList(sortIdx);
		
		% store event list
		save([dataPath eventPath eventFilename], 'fileName','eventList','eventTypes','raterID');
	end
	function cancelEvent(~, ~)
		% close event window
		delete(boxObj)
		boxObj = [];
		close(fEvent);
	end
	function updateEventListbox(eventIdx)
		nbEvents = numel(eventList);
		evStrings = cell(nbEvents, 1);
		for evIdx = 1:nbEvents
			% assemble event infos
			currEv = eventList{evIdx};
			evTime = currEv{2};
			evType = eventTypes{find(currEv{1}, 1)}(1:8);
			if(sum(currEv{1}) > 1)
				evType = strcat(evType, '/...');
			end
			evChan = chLabel{currEv{4}(1)};
			if(numel(currEv{4}) > 1)
				evChan = strcat(evChan, '/...');
			end
			
			evStrings{evIdx} = [num2str(evTime) ': ' evType ' (' evChan ')'];
			
			% edit event label appearance
			useHTML = [~any(selEvType(currEv{1})) isempty(intersect(currEv{4}, dispChan)) ~isempty(currEv{5})];
			if(useHTML(1) || useHTML(2))
				% grayed out if event type not selected or no event channel displayed
				evStrings{evIdx} = ['<font color="gray">' evStrings{evIdx} '</font>'];
			end
			if(useHTML(3))
				% bold and italic if commented
				evStrings{evIdx} = ['<b><i>' evStrings{evIdx} '</i></b>'];
			end
			if(any(useHTML))
				evStrings{evIdx} = ['<html>' evStrings{evIdx} '</html>'];
			end
		end

		set(lbxEvents, 'String',evStrings, 'Value',eventIdx);
	end
	function delEvent(~, ~)
		% delete currently selected event
		evIdx = get(lbxEvents, 'Value');
		if(evIdx)
			eventList(evIdx) = [];
			save([dataPath eventPath eventFilename], 'fileName','eventList','eventTypes','raterID');
			updateEventListbox(min(evIdx, numel(eventList)));
			plotEEGsignals();
		end
	end
	function gotoEvent(source, ~)
		% change current position to selected event
		evIdx = source.Value;
		if(~isempty(evIdx))
			evTime = eventList{evIdx}{2};

			winStartIdx = max(evTime - winLength/4, 0);
			set(sldTime, 'Value',winStartIdx);
		end
		
		plotSpectrumPosition();
		plotEEGsignals();
	end
	function plotEvents()
		evStarts = cellfun(@(x) x{2}, eventList);
		evEnds = cellfun(@(x) x{3}, eventList);
		dispEvTime = evStarts < (winStartIdx + winLength) & evEnds > winStartIdx;
		
		% get events to display
		dispEvType = cellfun(@(x) any(x{1}(selEvType)), eventList);
		dispEvIdx = dispEvTime & dispEvType;
		dispEv = eventList(dispEvIdx);
		
		% define colors of event boxes
		[~, ~, dispEvSel] = intersect(get(lbxEvents, 'Value'), find(dispEvIdx));
		evCols = repmat([0 1 0], numel(dispEv), 1);
		evCols(dispEvSel, :) = 1 - evCols(dispEvSel, :);
		
		for evIdx = 1:numel(dispEv)
			evStart = dispEv{evIdx}{2};
			evEnd = dispEv{evIdx}{3};
			
			evChan = dispEv{evIdx}{4};
			[~, ~, dispEvChIdx] = intersect(evChan, dispChan);
			if(tbRefBI.Value)
				dispEvChIdx = bsxfun(@eq, dispEvChIdx', biChCombs(:, 1)) | ...
								bsxfun(@eq, dispEvChIdx', biChCombs(:, 2));
				dispEvChIdx = find(any(dispEvChIdx, 2));
			end
			
			% draw colored boxes
			for chIdx = 1:numel(dispEvChIdx)
				minLen = 0.05;	% minimal box length
				xLen = evEnd - evStart;
				if(xLen < minLen)
					xStart = evStart - (minLen - xLen) / 2;
					xLen = minLen;
				else
					xStart = evStart;
				end
				
				yStart = dispEvChIdx(chIdx) * offset - offset/2;
				
				rectangle('Position',[xStart yStart xLen offset], ...
							'EdgeColor',[evCols(evIdx, :) 0.7], ...
							'FaceColor', [evCols(evIdx, :) 0.2], 'Parent',axEEG, ...
							'ButtonDownFcn',@modifyEvent);
				
				% add label on upper border(s)
				if(chIdx == 1 || dispEvChIdx(chIdx) > (dispEvChIdx(max(chIdx-1, 1))+1))
					text(max(xStart, winStartIdx+winLength/100), yStart-offset/3, ...
						strjoin(eventTypes(dispEv{evIdx}{1}), ' / '), ...
						'FontWeight','bold', 'Parent',axEEG, 'HitTest','off', 'Clipping','on');
				end
			end
		end
	end
	function modifyEvent(source, event)
		if(event.Button ~= 1)	% only continue on left click
			return
		end
		if(isvalid(fEvent))
			close(fEvent);		% close event window if opened already
		end
		
		% get idx of selected event based on time
		selEvStart = source.Position(1);
		selEvEnd = selEvStart + source.Position(3);
		evStarts = cellfun(@(x) x{2}, eventList);
		evEnds = cellfun(@(x) x{3}, eventList);
		evIdx = find((evStarts - selEvStart) <= 0.1 & abs(selEvEnd - evEnds) <= 0.1);
		
		% if ambiguous use channels to select
		if(sum(evIdx) > 1)
			selCh = ceil(source.Position(2) / offset);
			if(tbRefBI.Value)
				selCh = biChCombs(selCh, [2 1]);
			end
		
			selEvIdx = cellfun(@(x) ~isempty(intersect(dispChan(selCh), x{4})), eventList(evIdx));
			evIdx = evIdx(selEvIdx);
		end
		
		% if overlapping events, use closest centered
		if(sum(evIdx) > 1)
			evMids = evStarts(evIdx) + (evEnds(evIdx) - evStarts(evIdx)) / 2;
			evDist = abs(evMids - event.IntersectionPoint(1));
			[~, minDistIdx] = min(evDist);
			evIdx = evIdx(minDistIdx);
		end
		
		set(lbxEvents, 'Value',evIdx);
		plotEEGsignals();
		eventWin(evIdx);
	end
	
	
	%%% keyboard controls
	function keyDownListener(~, event, ~)
		switch(event.Key)
			case 'leftarrow'
				% go backward by tenth of window length
				currWinIdx = get(sldTime, 'Value');
				winStartIdx = max(currWinIdx - winLength/10, sldMin);
				set(sldTime, 'Value',winStartIdx);
			case 'rightarrow'
				% go forward by tenth of window length
				currWinIdx = get(sldTime, 'Value');
				winStartIdx = min(currWinIdx + winLength/10, sldMax);
				set(sldTime, 'Value',winStartIdx);
			case 'pageup'
				% go backward by window length
				currWinIdx = get(sldTime, 'Value');
				winStartIdx = max(currWinIdx - winLength, sldMin);
				set(sldTime, 'Value',winStartIdx);
			case 'pagedown'
				% go forward by window length
				currWinIdx = get(sldTime, 'Value');
				winStartIdx = min(currWinIdx + winLength, sldMax);
				set(sldTime, 'Value',winStartIdx);
			case 'home'
				% go to beginning of recording
				winStartIdx = 0;
				set(sldTime, 'Value',winStartIdx);
			case 'end'
				% go to end of recording
				winStartIdx = sldMax - winLength;
				set(sldTime, 'Value',winStartIdx);
			case 'uparrow'
				% increase amplitudes by factor 1.5
				if(~isempty(boxCoord))
					boxCoord(2) = boxCoord(2) / 1.5;
				end
				offset = offset / 1.5;
			case 'downarrow'
				% decrease amplitudes by factor 1.5
				if(~isempty(boxCoord))
					boxCoord(2) = boxCoord(2) * 1.5;
				end
				offset = offset * 1.5;
			case 'r'
				% zoom in by factor 2
				zoomHelpFun(max(winLength / 2, 0.1));
			case 't'
				% zoom out by factor 2
				zoomHelpFun(min(winLength * 2, sldMax));
			case 'e'
				% reset window to length 10s
				zoomHelpFun(10);
			case 'f'
				% switch filter on/off
				set(tbFilt, 'Value',1 - get(tbFilt, 'Value'));
				switchFilter();
			case 'c'
				% set record reference
				set(tbRefRec, 'Value',1);
				changeEEGsigDisp();
			case 'v'
				% set common average reference
				set(tbRefCA, 'Value',1);
				changeEEGsigDisp();
			case 'b'
				% set bipolar reference
				set(tbRefBI, 'Value',1);
				changeEEGsigDisp();
		end
		
		plotSpectrumPosition();
		plotEEGsignals();
	end
	

	%%% helper functions
	function zoomHelpFun(newWinLen)
		% set window start & length after zoom with constant center
		winMid = winStartIdx + winLength/2;
		winLength = newWinLen;
		winStartIdx = max(winMid - winLength/2, sldMin);
		winStartIdx = min(winStartIdx, sldMax - winLength);
		set(edtwinL, 'String',winLength);
		set(sldTime, 'SliderStep',[winLength/10 winLength] * fs / nbSamp);
		set(sldTime, 'Value',winStartIdx);
	end
	function setEventfileNameTxt(maxLen)
		if(numel(eventFilename) > maxLen)
			efnShort = ['...' eventFilename(max(1, numel(eventFilename)-maxLen+4):end)];
		else
			efnShort = eventFilename;
		end
		set(txtEvents, 'String',efnShort, 'Tooltip',eventFilename);
	end
	
end

