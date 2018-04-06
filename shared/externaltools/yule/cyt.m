% ------------------------
% SightOf is tool to analyze and visualize single cell data. Typically in
% situations where there are more cells than features (flow (7-10 features,
% mass-cytometry (30-40 features), microscopy, and more.
%
% SightOf Provides functionlity such as gating, clustering, 
% dimentionality reduction and comparing expression accross samples.
%
% ------
% 
% feel free to tweak the code.
% SightOf is 'event triggered'. Every GUI element you click triggers a
% callback function. Generally, functions you will be interested in are
% 'plotScatter','plot_histograms', 'runTSNE', 'calcModularity'. 
%
% Only two major global variables:
%
% 'sessionData', MXN numeric matrix. All the data loaded to SightOf. 
% Each fcs is concatinated at the bottom of the matrix. each column is a
% channel.
%
% 'gates', KX4 cell-matrix. Each row is a gate (K gates). 
% gates{i, 1} = visual name of gate that appears in the GUI
% gates{i, 2} = indices of the gate's data in the sessionData matrix
% gates{i, 3} = a horazontal cell array of channel names
% gates{i, 4} = filename if this gate was loaded from a file or empty
% 
% vairables are saved using "put('<varname>', var)"
% vairables are retrieved by "retr('<varname>')"
%
% Many place in the code use 'gateContext' which is all indices currently 
% selected. This is calculated whenever a user selects or deselects a
% gate from the gate list. 
%
% All UI components can be accessed from a handles global variable which
% is retrieved using the call 'gethand'.
%
% For example, a function that plots the selected channels would look like
% this:
%
% handles = gethand;
% sessionData = retr('SessionData');
% gateContext = retr('gateContext');
% selectedChannels = get(handles.lstChannels, 'Value');
% 
% plot(sessionData(gateContext, selectedChannels));
%
% ------ FAQ
%
% After a clustering or dim-reduction action the results are tacked on to
% the gates as extra channels.
% 
% Why do I have 'gate_source' Channels?
% When merging or subsampling from a number of gates to a single gate, a
% numeric 'gate_source' channel is added that you can later use to color
% by or to re-separate the data. 
%
% Why do I have 'cyt_placeholder' Channels?
% Channels are added as extra dimentions to the matrix as to not overwrite 
% another gate's channels. In some cases a 'cyt_placeholder' channels are 
% added to gates to achieve this. 
% This can be explained in the following example: 
% You gate all the high CD4 to a gate 'CD4+' and run tSNE. Now all the
% CD4 high cells have two extra channels. Then for some reason you want to
% tSNE all your 'T cells'. some of your T cells already have tSNE channels 
% so two placeholder channels are added to your 'T cells' gate before 
% running tSNE.
%
% Michelle Tadmor, Columbia University, 2012-2013
% ------------------------

function varargout = cyt(varargin)
% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @cyt_OpeningFcn, ...
                   'gui_OutputFcn',  @cyt_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT
end
% --- Executes just before cyt is made visible.
function cyt_OpeningFcn(hObject, ~, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to cyt (see VARARGIN)
    
    % initialize variable management
    handles.output = hObject;
    guidata(hObject, handles);
    handles=guihandles(hObject);
    setappdata(0,'hgui',gcf);


    % initialize visualizations per selected channels options.
    case0PlotTypes = {...	% plot options for all channels
        'Plot along CCT',        @plot_along_time};
%           'Histograms', 	@(by_gate)plot_histograms(1)};
%           'Histograms',             @(by_gate)plot_histograms(0)

    case1PlotTypes = {...	% plot options for a single selected channel
%           'Histograms',	@(by_gate)plot_histograms(1);
%           'Histograms',             @(by_gate)plot_histograms(0);
          'Plot along CCT',        @plot_along_time};
%           'Plot sample clusters',   @plot_sample_clusters;
%           'Plot meta clusters',     @plot_meta_clusters};

          
    case2PlotTypes = {... 	% plot options for 2 selected channels
          'Scatter',                @plotScatter;
%           'Density',                @plotDensity;
%           'Histograms',	@(by_gate)plot_histograms(1);
%           'Histograms',             @(by_gate)plot_histograms(0);
          'Plot along CCT',        @plot_along_time};
%           'Plot sample clusters',   @plot_sample_clusters;
%           'Plot meta clusters',     @plot_meta_clusters};
          
	case3PlotTypes = {...	% plot options for 3 selected channels
          'Scatter',                @plotScatter;
%           'Histograms',	@(by_gate)plot_histograms(1);
%           'Histograms',             @(by_gate)plot_histograms(0);
          'Plot along CCT',        @plot_along_time};
%           'Plot sample clusters',   @plot_sample_clusters;
%           'Plot meta clusters',     @plot_meta_clusters};
      
	case4PlusPlotTypes = {...	% plot options for 4 or more selected channels
%           'Histograms',	@(by_gate)plot_histograms(1);
%           'Histograms',             @(by_gate)plot_histograms(0);
          'Plot along CCT',        @plot_along_time};
%           'Plot sample clusters',   @plot_sample_clusters;
%           'Plot meta clusters',     @plot_meta_clusters};

    plotTypes = {case0PlotTypes, case1PlotTypes,...
        case2PlotTypes, case3PlotTypes, case4PlusPlotTypes};
    put('plotTypes', plotTypes);
    put('lastPlotTypeTable', ones(1, numel(plotTypes)));
    put('lastTimeChannel', 0);

    put('diff', 0);
    
    put('gates_listener',   @refreshGates)
      
    % if not found in settings - create default
    regexps = {'Cycler Features' '^Intensity_ResBlue_Nuclei_1_IntegratedIntensity$|^Intensity_ResFarRed_Nuclei_1_IntegratedIntensity$|^Texture_5_ResFarRed_Nuclei_12_InfoMeas1$|^Texture_5_ResBlue_Nuclei_4_Variance$|^Texture_5_ResBlue_Nuclei_7_SumVariance$';
               'Phase Indicators' '^G1$|^S$|^G2$';
               'All' '.*'};
    set(handles.lstRegexps, 'String', regexps(:, 1));
	put('regexps', regexps);

    % Update handles structure
    guidata(hObject, handles);
    
    % set defaults
    set(0,'DefaultTextInterpreter','None');
    
end

function put(name, what)
    hgui=getappdata(0,'hgui');
    setappdata(hgui, name, what);
    
    listener = retr([name '_listener']);
    if ~isempty(listener)
        listener();
    end
end

% -- returns saved variable or empty matrix if the variable is not found.
function var = retr(name)
    hgui=getappdata(0,'hgui');  
    var=getappdata(hgui, name);
end

function handles=gethand
    hgui=getappdata(0,'hgui');
    handles=guihandles(hgui);
end

% --- Outputs from this function are returned to the command line.
function varargout = cyt_OutputFcn(~, ~, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles;
end

% --------------------------------------------------------------------
function FileMenu_Callback(~, ~, ~)
% hObject    handle to FileMenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
end

function save_session
sessionData = retr('sessionData');
gates = retr('gates');
regexps = retr('regexps');
if (isempty(sessionData)) 
    uiwait(msgbox('The session is empty.','Nothing to save','modal'));
    return;
end

[filename,pathname,~] = uiputfile('*.mat','Save Session');

if isequal(filename,0) || isequal(pathname,0)
    return;
end

save([pathname filename], 'sessionData', 'gates', 'regexps'); 
%save([pathname filename], 'sessionData', 'gates', 'regexps', '-v7.3'); %alternative for data > 2GB

end

function load_session
handles = gethand;

% test if we are currently in session
sessionData = retr('sessionData');
if (~isempty(sessionData)) 
    selection = questdlg('This will remove your current session and load your selected one. Changes to your current session will not be saved. Are you sure you want to go ahead with this?' ,...
                         'Load SightOff Session',...
                         'Yes','No','No');
    if strcmp(selection,'No')
        return;
    end
end
clear sessionData;

% request session file
[filename, pathname, ~] = uigetfile('*.mat', 'Load SightOff Session');
if isequal(filename,0) || isequal(pathname,0)
    return;
end

% load file
load([pathname filename]);

% test file validity
if (exist('session_data', 'var')) %backwards compatibility
    sessionData = session_data;
    clear session_data;
end
if (~exist('sessionData', 'var') || ~exist('gates', 'var'))
    uiwait(msgbox('Unable to load needed data: ''SessionData'' or ''gates''. Check The file was saved by SightOf and contact the SightOf dev team with the file attatched for further investigation.','Load failed.','modal'));    
    return;
end

% set data
put('sessionData', sessionData);
put('gates', gates);

if exist('regexps', 'var')
    old_regexps = retr('regexps');
    for i = 1:size(regexps, 1)
        idx = strcmp(regexps{i, 1}, old_regexps(:,1));
        if ~any(idx)
            old_regexps{end+1, 1} = regexps{i, 1};
            old_regexps{end, 2} = regexps{i, 2};
            put('regexps', old_regexps);
            set(handles.lstRegexps, 'String', old_regexps(:, 1));
        end
    end
end

% make sure axis popups are visible and filled out
set(handles.plotChannels, 'Visible','on','Enable','on');
set(handles.pupPlotType, 'Visible','on','Enable','on');
 
selectedGate = get(handles.lstGates, 'Value');
if (selectedGate > size(gates, 1))
    set(handles.lstGates, 'Value', 1);
end
lstGates_Callback;

end

% --------------------------------------------------------------------
function OpenMenuItem_Callback(~, ~, ~)
    handles = gethand;
    
    currentfolder = retr('currentfolder');
    if (isempty(currentfolder)) 
        currentfolder = '';
    end

    files = uipickfiles('num',1,'out','cell', 'FilterSpec', [currentfolder '*.mat']);
    if isequal(files,0) ~=0
        return 
    end
    
    [path, ~, ext] = fileparts(files{1});
    put('currentfolder', [path filesep]);
    
    if (strcmp(ext, '.txt')) 
        % assume all files were txt based and we're loading from SCRATCH
        hWaitbar = waitbar(0,'reading txt files ...');
        tic
        fcsdats = cellfun(@dlmread, files, 'UniformOutput', false);
        disp(sprintf('Files loaded: %gs',toc));

        sessionData = zeros(0, size(fcsdats{1}, 2));
        gates = cell(numel(fcsdats),4);
        waitbar(0.5, hWaitbar, 'Adding data to session ...')
        for i=1:numel(fcsdats)
            
            [~, fcsname, ~] = fileparts(files{i}); 
            
            %-- add data to giant matrix
            currInd = size(sessionData, 1);
            sessionData(currInd+1:currInd+size(fcsdats{i},1), 1:size(fcsdats{i},2)) = fcsdats{i}(:, :);
            
            gates{i, 1} = char(fcsname);
            gates{i, 2} = currInd+1:currInd+size(fcsdats{i},1);
            gates{i, 3} = strcat({'channel '},int2str((1:size(fcsdats{i},2)).'))';
            gates{i, 4} = files{i}; % opt cell column to hold filename
        end
    elseif (strcmp(ext, '.h5')) % for prisca: assume only 1 file, empty session data
        [~, fcsname, ~] = fileparts(files{1}); 
        tic;
        hWaitbar = waitbar(0,'reading h5 file ...');
        info = h5info(files{1});
        sessionData = h5read(files{1}, [info.Name info.Datasets.Name]);
        disp(sprintf('File loaded: %gs',toc));
        
        gates = cell(1,4);
        gates{1, 1} = char(fcsname);
        gates{1, 2} = 1:size(sessionData,1);
        gates{1, 3} = strcat({'channel '},int2str((1:size(sessionData,2)).'))';
        gates{1, 4} = files{1}; % opt cell column to hold filename
    elseif (strcmp(ext, '.csv'))

        hWaitbar = waitbar(0,'reading csv files...');
        tic
        
        
        %Read in header
        fids = cellfun(@(f) fopen(f, 'r'), files, 'UniformOutput', false);  %opening files
        csvheaders = cellfun(@(f) fgetl(f), fids, 'UniformOutput', false);   %readinf first line in file
        cellfun(@(f) fclose(f), fids, 'UniformOutput', false);  %closing files
        
        %Convert  header to cell array
        csvheaders = cellfun(@(h) regexp(h, '([^,]*)', 'tokens'), csvheaders, 'UniformOutput', false);
        csvheaders = cellfun(@(h) cat(2, h{:}), csvheaders, 'UniformOutput', false);
        

        %Read in data    
        csvdats = cellfun(@(fname) csvread(fname, 1,0), files, 'UniformOutput', false);
        
        %Add data to session
        disp(sprintf('Files loaded: %gs',toc));

        sessionData = zeros(0, size(csvdats{1}, 2));
        gates = cell(numel(csvdats),4);
        waitbar(0.5, hWaitbar, 'Adding data to session ...')
        for i=1:numel(csvdats)
            
            [~, csvname, ~] = fileparts(files{i}); 
            
            %-- add data to giant matrix
            currInd = size(sessionData, 1);
            sessionData(currInd+1:currInd+size(csvdats{i},1), 1:size(csvdats{i},2)) = csvdats{i}(:, :);
            
            gates{i, 1} = char(csvname);
            gates{i, 2} = currInd+1:currInd+size(csvdats{i},1);
            gates{i, 3} = csvheaders{i};
%             gates{i, 3} = strcat({'channel '},int2str((1:size(csvdats{i},2)).'))';
            gates{i, 4} = files{i}; % opt cell column to hold filename
        end
        
    elseif (strcmp(ext, '.mat')) % for Gabriele: assume only 1 file, empty session data
           [~, fcsname, ~] = fileparts(files{1}); 
           tic;
           hWaitbar = waitbar(0,'reading mat file ...');
           load(files{1}); % assume now we have data and featurenames
           sessionData = data;

           disp(sprintf('File loaded: %gs',toc));

           gates = cell(1,4);
           gates{1, 1} = char(fcsname);
           gates{1, 2} = 1:size(sessionData,1);
           gates{1, 3} = featurenames;
           gates{1, 4} = files{1}; % opt cell column to hold filename
    else% assume files are fcs format (the only officialy supported format
        
        tic
        [fcsdats fcshdrs]=cellfun(@fca_readfcs, files, 'UniformOutput', false);
        %cdatas = cellfun(@cytof_data, files, 'UniformOutput', false );

        disp(sprintf('Files loaded: %gs',toc));

        tic
        % find out how many 'channels' to allocated
        y = 0;
        nfcs = size(fcsdats, 2);
        hWaitbar = waitbar(0,'Allocating space for session data ...');
        for i=1:nfcs
            waitbar(i/nfcs, hWaitbar);
            y = max([y size(fcsdats{i}, 2)]);
        end

        % read all data to one huge matrix
        % and defined gates according to each filename
        sessionData  = retr('sessionData');
        if (isempty(sessionData)) 
            sessionData = zeros(0, y);
            gates = cell(nfcs,4);
            last_gate_ind = 0;
        else 
            gates = retr('gates');
            last_gate_ind = size(gates, 1);

            % if we're adding gates that have extra channels. like after the
            % user has ran tSNE or something like that
            if (size(sessionData, 2)< y) 
                sessionData(:, end+1:y) = zeros(size(sessionData,1), y - size(sessionData,2));
            end
        end
        disp(sprintf('Allocated space for data: %gs',toc));

        tic
        waitbar(0, hWaitbar, 'Adding data to session ...')
        for i=1:nfcs

            %-- add data to giant matrix
            currInd = size(sessionData, 1);
            sessionData(currInd+1:currInd+size(fcsdats{i},1), 1:size(fcsdats{i},2)) = fcsdats{i}(:, :);

            %-- save files as gates
            [~, fcsname, ~] = fileparts(files{i}); 
            gates{last_gate_ind+i, 1} = char(fcsname);
            gates{last_gate_ind+i, 2} = currInd+1:currInd+size(fcsdats{i},1);
%            gates{last_gate_ind+i, 3} = cdatas{i}.channel_name_map;        
             gates{last_gate_ind+i, 3} = get_channelnames_from_header(fcshdrs{i});        
            gates{last_gate_ind+i, 4} = files{i}; % opt cell column to hold filename

            waitbar(i-1/nfcs, hWaitbar, sprintf('Adding %s data to session  ...', gates{last_gate_ind+i, 1}));
        end
        disp(sprintf('Read data into session: %gs',toc));

    end
	waitbar(1, hWaitbar, 'Saving ...');

    % save the huge matrix, gates and use default empty gate
    put('sessionData', sessionData);
    put('gates', gates);
    
    set(handles.lstGates,'Value',[]); % Add this line so that the list can be changed
    set(handles.lstGates,'String',gates(:, 1));
    set(handles.lstIntGates,'Value',[]); % Add this line so that the list can be changed
    set(handles.lstIntGates,'String',gates(:, 1));
    
    set(handles.lstChannels,'Value',[]); % Add this line so that the list can be changed
    
    % make sure axis popups are visible and filled out
    set(handles.plotChannels, 'Visible','on','Enable','on');
    set(handles.pupPlotType, 'Visible','on','Enable','on');
    
    set(handles.lstGates,'Value',1);
    set(handles.lstChannels,'Value',[]);
    
    lstGates_Callback;
    lstChannels_Callback;
    
    close(hWaitbar);
end

function channel_names=get_channelnames_from_header(fcshdr)
    channel_names = {fcshdr.par.name2};
	if isempty(channel_names{1})
        channel_names = {fcshdr.par.name};
	end
end

% --------------------------------------------------------------------
function PrintMenuItem_Callback
    [filename, pathname, ~] = uiputfile({'*.png;*.pdf', 'Image files (.png, .pdf)';...
                                         '*.png', 'Portable Network Graphics (.png)';...
                                         '*.pdf', 'Color PDF file format (.pdf)'},...
                                        'Save Image');

    if isequal(filename,0) || isequal(pathname,0)
    return;
    end
    
    ha = gca;
    f_new = figure;
    copyobj(ha, f_new);
    
    % save new figure
    if (endswith(filename, '.pdf'))
        screen2eps(f_new, [pathname filename]);
    else
        screen2png(f_new, [pathname filename]);
    end
end
% --------------------------------------------------------------------
function CloseMenuItem_Callback(~, ~, handles)
% hObject    handle to CloseMenuItem (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
selection = questdlg(['Close ' get(handles.figure1,'Name') '?'],...
                     ['You will lose any unsaved work.\nClose ' get(handles.figure1,'Name') '...'],...
                     'Yes','No','Yes');
if strcmp(selection,'No')
    return;
end

delete(handles.figure1)
end

% --- Executes on selection change in lstChannels.
function lstChannels_Callback(~, ~, ~)
    handles = gethand;

    % -ignore- 
    % Workaround for bug where context menus don't show on mac.
    ctm = get(handles.lstChannels, 'UIContextMenu'); 
    % -ignore- 

    % -- save newly selected channels
    selectedChannels  = get(handles.lstChannels,'Value');
    nSelectedChannels = numel(selectedChannels);

    % -- update plot options according the # of selected channels.
    plotTypes         = retr('plotTypes');
    lastPlotTypeTable = retr('lastPlotTypeTable');
    nSelectedType     = min(nSelectedChannels+1, numel(plotTypes));
    plotBynSelected   = plotTypes{nSelectedType};
    set(handles.pupPlotType, 'String', plotBynSelected(:,1), 'Value', lastPlotTypeTable(nSelectedType));

    % -- clear regexps textfield
    set(handles.txtRegexps, 'String', '');
end

% --- Executes on selection change in lstGates.
function lstGates_Callback(~, ~, ~)
    handles = gethand;

    % -ignore- matlab right click menu not showing bug workaround
    ctm = get(handles.lstGates, 'UIContextMenu'); 
    
    gates = retr('gates');
    
    if (size(gates, 1) == 0)
        return;
    end
    

    selected_gates = get(handles.lstGates, 'Value');
    selected_int_gates = get(handles.lstIntGates, 'Value');
    
    % if no gate is selected, then the context is all the data available 
    % in the session.
    if isempty(selected_gates)
        selected_gates = 1:size(gates, 1);
    end
    
    [gate_indices, channel_names] = getSelectedIndices(selected_gates);

    % filter indices if user selected to intersect gates
    if get(handles.btnIntersect, 'Value')
        [gate_int_indices channel_int_names] = getSelectedIndices(selected_int_gates);
        if (~isempty(gate_int_indices))
            gate_indices = intersect(gate_int_indices, gate_indices);
            if (numel(channel_int_names) > numel(channel_names))
                channel_names = channel_int_names;
            end
        end    
    end    
            
    % save new channel names
    put('channelNames', channel_names);
    
    % attempt to preserve selection
    selected_channels = get(handles.lstChannels, 'Value');
    selected_channels(find(selected_channels > numel(channel_names))) = [];
    set(handles.lstChannels, 'Value', selected_channels);

    % set gui components with new channels names
    set(handles.lstChannels, 'String', channel_names);
    setChannelNamesToAxis(handles.lstClusterChannels, channel_names);
    setChannelNamesToAxis(handles.lstMetaClusterChannels, channel_names);
    setChannelNamesToAxis(handles.lstBasisOfMetaChannel, channel_names);
    setChannelNamesToAxis(handles.lstTsneColorBy, channel_names);
    setChannelNamesToAxis(handles.lstFunctionalXAxis, channel_names);
    setChannelNamesToAxis(handles.pupXAxis, channel_names);
    setChannelNamesToAxis(handles.pupYAxis, channel_names);
    setChannelNamesToAxis(handles.pupDensXAxis, channel_names);
    setChannelNamesToAxis(handles.pupDensYAxis, channel_names);
    setChannelNamesToAxis(handles.pupZAxis, [' ' channel_names(:)']);
    setChannelNamesToAxis(handles.pupColorBy,['Gates' channel_names(:)']);

    
    % show\hide option to plot difference of channels (TODO move this somewhere
    % else!!)
    if (numel(selected_gates)== 2)
        %You can add in the future 
        set(handles.btnDiff, 'Visible', 'off');
%         set(handles.btnDiff, 'Visible', 'on');
%         setChannelNamesToAxis(handles.lstKNNPresentationSpace, gates{selected_gates(1), 3});
%         setChannelNamesToAxis(handles.lstKNNSpace, channel_names);
    else
        set(handles.btnDiff, 'Visible', 'off');
    end
    
    % save new gate context
    put('gateContext', gate_indices);
    
    setStatus(sprintf('%g data points in the gates selected', numel(gate_indices)));
end

function setChannelNamesToAxis(pupObject, channel_names)
    set(pupObject, 'String', channel_names, 'Value', min(numel(channel_names), get(pupObject, 'Value')));
end

function [indices channels] = getSelectedIndices(selected_gates)
    gates = retr('gates');
    
    % extract specific gate or merge multiple gates according to selection
    if (numel(selected_gates) == 1)
        indices = gates{selected_gates, 2};
        channels =  gates{selected_gates, 3};
    else 
        indices = [];
        if (~isempty(selected_gates))
            channels = gates{selected_gates(1),3};
        else
            channels = [];
        end
        
        % --- for simplicity we assume same channels for all gates  except for
        % trailing channels. so changes in size are the only changes in channels. 
        
        % loop thorugh each selected gate. we'll collect (union) the data
        % and 'intersect' the channels as some gates may have more or
        % different channels appended.
        for i=selected_gates
            
            indices = union(gates{i, 2}, indices);            
            
            if (size(channels,2) > size(gates{i,3}, 2)) 
                
                % shorten the channel names in use
                channels = channels(1:size(gates{i,3}, 2));
            end
        end
    end
end

function intIndices=getIntIndices
    handles = gethand;
    gates = retr('gates');
    intIndices = [];
    if (get(handles.btnIntersect, 'Value'))
        selIntGates = get(handles.lstIntGates, 'Value');
        for i=selIntGates
            intIndices = union(intIndices, gates{i, 2});
        end
    end
end
% --- Executes on button press in btnAddFCS.
function btnAddFCS_Callback(~, ~, ~)
% hObject    handle to btnAddFCS (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
end

% --- Executes on button press in plotChannels.
function plotChannels_Callback(~, ~, handles)
    
    % disable gating functionality
    set(handles.btnGatePoly, 'Enable', 'off');
    set(handles.btnGateRect, 'Enable', 'off');
    
    % find selected plot type request
    selectedGates = get(handles.lstGates,'Value');
    selectedChannels = get(handles.lstChannels,'Value');
    nSelectedChannels = numel(selectedChannels);
    idxSelectedPlotType = get(handles.pupPlotType, 'Value');
    plotTypes        = retr('plotTypes');
    selectedPlotType = retr('lastPlotTypeTable');
    
    nSelectedChannelsType = min(nSelectedChannels+1, numel(plotTypes));
    plotBynSelected  = plotTypes{nSelectedChannelsType};
    
    % Open in a New Window if Command/ctrl is pressed
	current_figure = gcf;
    if (isCtrlPressed)
        h = figure('Color',[1 1 1]);
    else
        % hide panels
        hidePlotControls;
    end
    
    % Save selected gates, channels, and plot type
    put('inMainPlot', true);

    % Save selected gates, channels, and plot type
    put('selectedGates', selectedGates);

    put('selectedChannels', selectedChannels);

    selectedPlotType(nSelectedChannelsType) = idxSelectedPlotType;
    put('lastPlotTypeTable', selectedPlotType);

    
    % Invoke the selected plot function
    plotBynSelected{idxSelectedPlotType, 2}();

    % Save selected gates, channels, and plot type
    put('inMainPlot', false);
 
    % Return to base figure
    set(0,'CurrentFigure', current_figure);
end

function isCtrlPressed=isCtrlPressed
    modifiers = get(gcf,'currentModifier'); 
    isCtrlPressed = ~isempty(find(ismember({'command'; 'control'},modifiers)));
end

function hidePlotControls
    handles = gethand;


    set(handles.pnlScatterControls, 'Visible', 'off');
    set(handles.pnlDensityControls, 'Visible', 'off');
    set(handles.pnlHistogramControls, 'Visible', 'off');
    set(handles.pnlWonderlustControls, 'Visible', 'off');
    set(handles.pnlClusterControls, 'Visible', 'off');
    set(handles.pnlSingleClusterControls, 'Visible', 'off');
    set(handles.pnlMetaClusterControls, 'Visible', 'off');
    set(handles.pnlSingleMetaControls, 'Visible', 'off');
    set(handles.pnlMetaTsneColor, 'visible', 'off');
    
end

function plot_histograms(by_gates) 
    % ---- init variables and clear prepare plotting planel ----

    handles = gethand;

    sessionData	= retr('sessionData');
    gateContext	= retr('gateContext');
    intIndices  = getIntIndices;
	gates       = retr('gates');
    channelNames = retr('channelNames');

    selectedChannels = get(handles.lstChannels,'Value');
    nSelChannels     = numel(selectedChannels);

    selectedGates = retr('selectedGates');
    displayChannels = selectedChannels;
    nSelGates     = numel(selectedGates);
    set(handles.lstHistGroupA, 'Value', selectedGates);
    
    % clear the figure panel
    cla;
    colormap jet;
    legend('off');
    colorbar off;
    
    if (isempty(gateContext)) return; end
    
    subplot(1,1,1,'Parent',handles.pnlPlotFigure);

    if (nSelChannels <= 0) 
        return;
    end
    
    if by_gates && numel(selectedGates) > 1
    end

    % fix axis labels and popups
    if ~isCtrlPressed
        hidePlotControls
        set(handles.pnlHistogramControls, 'Visible','on');
        if ~isempty(get(handles.lstHistGroupB, 'Value'))
            set(handles.pnlHistogramTop, 'Visible','on');
        else
            set(handles.pnlHistogramTop, 'Visible','off');
        end
   end
    
    % ---- begin plotting ----
    nrows = round(sqrt(nSelChannels+1));

    if nrows^2 < (nSelChannels+1)
        ncols = nrows+1;
    else
        ncols = nrows;
    end
    
    if nSelChannels == 2
        nrows = 2;
        ncols = 3;
    end

    % set colors and reset num of subplots
    ColorMtx = distinguishable_colors(max(2, nSelGates));
    subplot(1,1,1,'Parent',handles.pnlPlotFigure);
    dplot(1:100, 'colors', ColorMtx);
    legend_space = 0;

    % plot legend
    if (by_gates)

        legend_space = 1; % flag to present a legend
        
        subpop_selection = get(handles.lstHistGroupB, 'Value');
        if ~isempty(subpop_selection)
            
            % -- sort 'selected gates' by difference           
            refpop_selection = get(handles.lstHistGroupA, 'Value');
            refInd = getSelectedIndices(refpop_selection);
            subInd = getSelectedIndices(subpop_selection);
            ref_sub_intInd = intersect(refInd, subInd);
            if ~get(handles.btnHistIntLeft, 'Value')
                refInd = setdiff(refInd, ref_sub_intInd);
            end
            if ~get(handles.btnHistIntRight, 'Value')
                subInd = setdiff(subInd, ref_sub_intInd);                
            end
            if ~isempty(intIndices)
                refInd = intersect(refInd, intIndices);
                subInd = intersect(subInd, intIndices);
            end
            
            referenceData = sessionData(refInd, :);
            subpopData    = sessionData(subInd, :);
            
            if isempty(referenceData) || isempty(subpopData)
                uiwait(msgbox('One of the populations is contained in the other, please allow the intersection of the groups to be accounted by selecting the intersection arrow below, otherwise it is empty and its histogram is not shown.','Cannot compare populations.','error'));               
                return;
            end
            
            hwaitbar = waitbar(0,'comparing distributions ...');
            
            for j = 1:nSelChannels
                dists(j) = histDist( referenceData(:,selectedChannels(j)) , subpopData(:,selectedChannels(j)) );
                waitbar(j/nSelChannels,hwaitbar, 'comparing distributions ...');
            end
           
            [~,ix] = sort( abs(dists), 'descend');
            
            waitbar(1,hwaitbar, 'Done.');

            close(hwaitbar);

            top = str2num(get(handles.txtHistTop, 'String'));
            if (top == 0) 
                top = numel(selectedChannels);
            else
                nrows = round(sqrt(top+1));

                if nrows^2 < (top+1)
                    ncols = nrows+1;
                else
                    ncols = nrows;
                end

                if top == 2
                    nrows = 2;
                    ncols = 3;
                end
            end
            displayChannels = selectedChannels(ix(1:min(top, numel(selectedChannels))));
            
            hPlot = subplot(nrows, ncols, 1);
            for j = 1:2 % generate two dummy hist lines for future legend
                hLine = dplot(1:100, 'colors', ColorMtx);
                set(hLine,'Visible','off');
            end
            legend({'reference group', 'subpopulation group'}, 'Interpreter', 'none');
        else
            % generate as many dummy lines as are gates selected
            hPlot = subplot(nrows, ncols, 1);
           for j = selectedGates
                hLine = dplot(1:100, 'colors', ColorMtx);
                set(hLine,'Visible','off');
            end
            legend(remove_repeating_strings(gates(selectedGates, 1)), 'Interpreter', 'none');          
        end
 
        set(hPlot,'Visible','off'); % hide dummy plot and keep legend ;)
    end

    % plot histograms
    i = 1+legend_space;
    for channel = displayChannels
        subplot(nrows, ncols, i);
        i = i+1;
        if (by_gates)
            if ~isempty(subpop_selection) 
                dplot(referenceData(:, channel),'colors',ColorMtx);
                dplot(subpopData(:, channel),'colors',ColorMtx);
            else
                for j = selectedGates
                    currGate = gates{j, 2};
                    if ~isempty(intIndices)
                        currGate = intersect(intIndices, currGate);
                    end
                    dplot(sessionData(currGate, channel),'colors',ColorMtx);
                end
            end
        else
            dplot(sessionData(gateContext, channel),'colors',ColorMtx);
        end
        box on
        set(gca,'ytick',[]);
        set(gca,'yticklabel',{});
        if exist('dists', 'var')
            title(sprintf('%s\ndiff: %g', channelNames{channel}, dists(find(channel==selectedChannels))));
        else
            title(channelNames(channel));
        end
    end

   sprintf('Plotted %i points', length(gateContext));

end

function plot_along_time(time_channel)
    handles = gethand; 

    % clear the figure panel
    cla;
    colormap jet;
    legend('off');
    colorbar('delete');
    axis auto;
    
    % show controls
    if ~isCtrlPressed
        set(handles.pnlWonderlustControls, 'Visible', 'on');
    end

	hPlot = subplot(1,1,1,'Parent',handles.pnlPlotFigure);
    box on;

    session_data = retr('sessionData'); % all data
    gates        = retr('gates');
    gate_context = retr('gateContext'); % indices currently selected
    if isempty(gate_context)
        return;
    end
    
    selected_channels = get(handles.lstChannels, 'Value');
    channel_names     = retr('channelNames');
    selected_gates    = get(handles.lstGates, 'Value');
    gate_names        = gates(selected_gates, 1);
    
    time_channel      = retr('lastTimeChannel');
%     time_channel      = get(handles.lstFunctionalXAxis, 'Value');

    if time_channel==0
        chcct = cellstrfnd(channel_names, 'cct');
        if any(chcct)
            s = find(chcct);
            s = s(1);
        else
            initTimeSelection = numel(channel_names);
            [s,v] = listdlg('PromptString','Select a Cell cycler trajectory feature:',...
                    'SelectionMode','single',...
                    'InitialValue', initTimeSelection,...
                    'ListString',channel_names);
            if ~v
                return;
            end
        end
        time_channel = s;
        set(handles.lstFunctionalXAxis, 'Value', time_channel);
        put('lastTimeChannel', time_channel);
    else
        time_channel = get(handles.lstFunctionalXAxis, 'Value');
    end
    
    % ==== TODO === ask by gate or by channels
    by_gate = false;
    
    arrWonderlust = session_data(gate_context, time_channel);
    matData       = session_data(gate_context, selected_channels);
    
    smoothness_factor   = str2num(get(handles.txtWindowSize, 'String'));
    normalizeX          = get(handles.chkWanderlustNormX, 'Value');
    normalizeY          = get(handles.chkWanderlustNormY, 'Value');
    rankY               = get(handles.chkRankY, 'Value');
    svGolay             = get(handles.chkSGolay, 'Value');
    show_error          = get(handles.chkWanderlustError, 'Value');

    avg_types           = get(handles.pupAvgType, 'String');
    avg_type            = avg_types{get(handles.pupAvgType, 'Value')};
    
    if by_gate
        % generate gate source vector for grouping
        for j=selected_gates
            v(ismember(gate_context,gates{j, 2})) = j;
        end

        % set legend labels to gate names
        legend_labels = gate_names;
    else
        v = [];
        legend_labels = channel_names(selected_channels);
    end

    if normalizeX
        % display wonderlust results
        arrWonderlust = arrWonderlust-min(arrWonderlust);
        arrWonderlust = arrWonderlust./max(arrWonderlust);
    end
    
    plot_as_function(arrWonderlust, matData, ...
                    'num_locs', 100,...
                    'avg_type', avg_type,...
                    'show_error', show_error,...
                    'labels', channel_names(selected_channels),...
                    'normalize', normalizeY,...
                    'rank', rankY,...
                    'svGolay', svGolay,...
                    'smooth', smoothness_factor);
    xlabel(channel_names{time_channel}); 
    if (numel(selected_channels) ==1)
        ylabel(channel_names{selected_channels(1)});
%         [mi, ~, ~, ~] = compute_dremi([arrWonderlust matData], .9);
%         title(sprintf('MI: %g', mi));
    else
        title('');
    end
%     legend(legend_labels);

    if (by_gate)
        title(channel_names(selected_channels(1)), 'Interpreter', 'none');
    end
end

function plot_sample_clusters(cluster_channel)
    handles = gethand; 

    % clear the figure panel
    cla;
    colormap jet;
    legend('off');
    colorbar('delete');
    axis auto;
    
    % show controls
    if ~isCtrlPressed
        set(handles.pnlClusterControls, 'Visible', 'on');
    end

	hPlot = subplot(1,1,1,'Parent',handles.pnlPlotFigure);
    box on;

    session_data = retr('sessionData'); % all data
    gates        = retr('gates');
    gate_context = retr('gateContext'); % indices currently selected
    if isempty(gate_context)
        return;
    end
    
    selected_channels = get(handles.lstChannels, 'Value');
    channel_names     = retr('channelNames');
    selected_gates    = get(handles.lstGates, 'Value');
    gate_names        = gates(selected_gates, 1);
    
    initClusterSelection = size(session_data, 2);
    
    %find cluster channel
    cluster_channel      = get(handles.lstClusterChannels, 'Value');
    if isempty(cluster_channel) || ~isDiscrete(cluster_channel) || length(unique(session_data(gate_context, cluster_channel))) > 500, %if channels is empty or channel is not discrete or the number of clusters is too large
%         [s,v] = listdlg('PromptString','Select a cluster channel:',...
%                 'SelectionMode','single',...
%                 'InitialValue', initClusterSelection,...
%                 'ListString',channel_names);
        if ~v
            return;
        elseif ~isDiscrete(cluster_channel)
            return;
        elseif length(unique(cluster_channel)) > 500,
            return;
        end
        cluster_channel = s;
    end

    
    %find plot method, 1=Heat map, 2=Pie chart, 3=tSNE
    %plot_type      = get(handles.popSamplePlotOptions, 'Value');
    
    %finding plot type
    plot_types = get(handles.popSamplePlotOptions, 'String');
    plot_type = plot_types(get(handles.popSamplePlotOptions, 'Value'));

    %plot types are:
%     Heat map: Mean
%     Heat map: Distance
%     Pie chart
%     tSNE map
%     Single cluster

    
    
    if strcmp(plot_type, 'Heat map: Mean')
        
        show_by_cluster_channel = true;
        if (show_by_cluster_channel)
            num_clusters = length(unique(session_data(gate_context, cluster_channel)));

            %finding mean values of marker levels for each cluster
            marker_means = zeros(num_clusters, length(selected_channels));
            data = session_data(gate_context, selected_channels);

            %looping through clusters
            clusters_in_sample = unique(session_data(gate_context, cluster_channel));
            for i=1:length(clusters_in_sample)
                marker_means(i,:) = mean(data(session_data(gate_context, cluster_channel)==clusters_in_sample(i),:),1);
            end
            marker_means = mynormalize(marker_means, 100);

            %find percentage of cells belonging to each cluster
            cells_pr_cluster = countmember(unique(session_data(gate_context, cluster_channel)),session_data(gate_context,cluster_channel))/length(gate_context);
        else
            num_clusters = length(selected_gates);

            %finding mean values of marker levels for each cluster
            marker_means = zeros(num_clusters, length(selected_channels));
            data = session_data(gate_context, selected_channels);

            %looping through gates
            for i=1:numel(selected_gates)
                marker_means(i,:) = mean(session_data(gates{selected_gates(i),2},selected_channels),1);
            end
            marker_means = mynormalize(marker_means, 100);

            %find percentage of cells belonging to each cluster
            cells_pr_cluster = cellfun(@(v)length(v), gates(selected_gates, 2))/length(gate_context);
        end

        %find color scheme
        color = get(handles.popSampleHeatColor, 'Value');

        %finding basis for min and max for color scaling
        heat_min = min(min(marker_means));
        heat_max = max(max(marker_means));
        
        if color == 1,
            color_map = interpolate_colormap(redbluecmap, 64);
            
            %adjusting heatmap color scaling
            if abs(heat_min) > heat_max,
                heat_max = abs(heat_min);
            elseif heat_min < 0,
                heat_min = -heat_max;
            end
        elseif color == 2,
            color_map = color_map_creator('rg');
            
            %adjusting heatmap color scaling
            if abs(heat_min) > heat_max,
                heat_max = abs(heat_min);
            elseif heat_min < 0,
                heat_min = -heat_max;
            end
        elseif color == 3,
            color_map = interpolate_colormap(othercolor('YlOrBr9'), 64);
        end
        
        %plot heat map
        %subplot(10,1,1:9), imagesc(hl_L2_dist, [heat_min,heat_max]);

        subplot(10,1,1:9), imagesc(marker_means);
        %set(title(strcat('Sample ', mat2str(p), ': L2 + higher/lower than median')));
        set(gca, 'Position', [0.18,0.27,0.8,0.66]);
        set(gca, 'ytick', 1:num_clusters);
        if (show_by_cluster_channel)
            y_labs = strcat(strread(num2str(unique(session_data(gate_context, cluster_channel))'), '%s'),' (',strread(num2str(cells_pr_cluster'), '%s'),')');
        else
            y_labs = strcat(gates(selected_gates, 1),' (',strread(num2str(cells_pr_cluster'), '%s'),')');            
        end
        set(gca, 'Yticklabel', y_labs);
        set(gca, 'xtick', []);
        colormap(color_map);
        colorbar;
        xticklabel_rotate([1:length(selected_channels)],90,channel_names(selected_channels));
        
        
    elseif strcmp(plot_type,'Heat map: Distance')
        %calculate data for heat map
        
        %using SPR
        
        
        %using L2 distance
        SPR_distances = SPR_dist_heatmap(session_data(gate_context, selected_channels), session_data(gate_context, cluster_channel), 3);

        %find percentage of cells belonging to each cluster
        cells_pr_cluster = countmember(unique(session_data(gate_context, cluster_channel)),session_data(gate_context,cluster_channel))/length(gate_context);
        
        heat_min = min(min(SPR_distances));
        heat_max = max(max(SPR_distances));
        
        %scaling colors so that zero will always be white
        if abs(heat_min) > heat_max,
            heat_max = abs(heat_min);
        elseif heat_min < 0,
            heat_min = -heat_max;
        end
        
        %defining color map
        color = get(handles.popSampleHeatColor, 'Value');
        
        if color == 1,
            color_map = interpolate_colormap(redbluecmap, 64);
            
            %adjusting heatmap color scaling
            if abs(heat_min) > heat_max,
                heat_max = abs(heat_min);
            elseif heat_min < 0,
                heat_min = -heat_max;
            end
        elseif color == 2,
            color_map = color_map_creator('rg');
            
            %adjusting heatmap color scaling
            if abs(heat_min) > heat_max,
                heat_max = abs(heat_min);
            elseif heat_min < 0,
                heat_min = -heat_max;
            end
        elseif color == 3,
            color_map = interpolate_colormap(othercolor('YlOrBr9'), 64);
        end
        
        %plot heat map
        subplot(10,1,1:9), imagesc(SPR_distances, [heat_min,heat_max]);
        %set(title(strcat('Sample ', mat2str(p), ': L2 + higher/lower than median')));
        set(gca, 'Position', [0.18,0.27,0.8,0.66]);
        set(gca, 'ytick', 1:length(unique(session_data(gate_context, cluster_channel))));
        y_labs = strcat(strread(num2str(unique(session_data(gate_context, cluster_channel))'), '%s'),' (',strread(num2str(cells_pr_cluster'), '%s'),')');
        set(gca, 'Yticklabel', y_labs);
        set(gca, 'xtick', []);
        colormap(color_map);
        colorbar;
        xticklabel_rotate([1:length(selected_channels)],90,channel_names(selected_channels));
        



%     elseif strcmp(plot_type,'Heat map: Distance')
%         %calculate data for heat map
%         
%        
%         
%         %using L2 distance
%         hl_L2_dist = L2_dist_heatmap_data(session_data(gate_context, selected_channels), session_data(gate_context, cluster_channel));
%         disp(sprintf('L2 dist heatmap computed: %gs',toc));
% 
%         %find percentage of cells belonging to each cluster
%         cells_pr_cluster = countmember(unique(session_data(gate_context, cluster_channel)),session_data(gate_context,cluster_channel))/length(gate_context);
%         
%         heat_min = min(min(hl_L2_dist));
%         heat_max = max(max(hl_L2_dist));
% 
%         %scaling colors so that zero will always be white
%         if abs(heat_min) > heat_max,
%             heat_max = abs(heat_min);
%         elseif heat_min < 0,
%             heat_min = -heat_max;
%         end
%         
%         %defining color map
%         color = get(handles.popSampleHeatColor, 'Value');
%         
%         if color == 1,
%             color_map = interpolate_colormap(redbluecmap, 64);
%         elseif color == 2,
%             color_map = color_map_creator('rg');
%         end
%         
%         
%         %plot heat map
%         subplot(10,1,1:9), imagesc(hl_L2_dist, [heat_min,heat_max]);
%         %set(title(strcat('Sample ', mat2str(p), ': L2 + higher/lower than median')));
%         set(gca, 'Position', [0.18,0.27,0.8,0.66]);
%         set(gca, 'ytick', 1:length(unique(session_data(gate_context, cluster_channel))));
%         y_labs = strcat(strread(num2str(unique(session_data(gate_context, cluster_channel))'), '%s'),' (',strread(num2str(cells_pr_cluster'), '%s'),')');
%         set(gca, 'Yticklabel', y_labs);
%         set(gca, 'xtick', []);
%         colormap(color_map);
%         colorbar;
%         xticklabel_rotate([1:length(selected_channels)],90,channel_names(selected_channels));
%     
    %if plot type is pie chart
    elseif strcmp(plot_type,'Pie chart'),
        
        %find percentage of cells belonging to each cluster
        cells_pr_cluster = countmember(unique(session_data(gate_context, cluster_channel)),session_data(gate_context,cluster_channel))/length(gate_context);
        
        labels = strread(num2str(unique(session_data(gate_context, cluster_channel))'),'%s');
        pie(cells_pr_cluster, zeros(1,length(cells_pr_cluster)), labels);
        colormap(distinguishable_colors(length(unique(session_data(gate_context,cluster_channel)))));
    
    %if plot type is tSNE map
    elseif strcmp(plot_type,'tSNE map')
        
        if length(gate_context) > 10000,
            inds = randsample(gate_context,10000);    %subsample if gate context is larger than 10,000
        else
            inds = gate_context;
        end
        
        clusters = session_data(inds, cluster_channel);
        
        perplexity = [];
        if(size(clusters,1) < 100)  %100 is arbitrary
            perplexity = 2;
            tsne_out = fast_tsne(session_data(inds,selected_channels),[], perplexity); %running tSNE
        else
            tsne_out = fast_tsne(session_data(inds,selected_channels),110); %running tSNE
        end
        
        %plotting scatterplot
        colors = distinguishable_colors(length(unique(clusters)));
        colormap(colors);
        clim = [min(clusters),max(clusters)];
        myplotclr(tsne_out(:,1), tsne_out(:,2), clusters, clusters, 'o', colors, clim,0);
        colorbar;
        caxis(clim);    %changing labels on colorbar
    
    %if plot type is single cluster
    elseif strcmp(plot_type, 'Single cluster')
        unique_clusters = unique(session_data(:,cluster_channel));
        
        if ismember(0,unique_clusters),
            fprintf('Cluster number 0 not allowed and will not be shown as an option\n')
            unique_clusters(unique_clusters == 0) = [];
        end
        
        set(handles.lstSingleCluster, 'String', unique_clusters);
        
        set(handles.txtSingleClusterDensTop, 'Value', 9);
        
        %channel_density_options = num2cell(1:length(selected_channels));
        %set(handles.popSingleCLusterDensityNumber, 'Value', channel_density_options);
        
        if ~isCtrlPressed
            set(handles.pnlSingleClusterControls, 'Visible', 'on');
            %set(handles.popSingleCLusterDensityNumber, 'Visible', 'on');

        end
        
        current_cluster = unique_clusters(get(handles.lstSingleCluster, 'Value'));
        plot_single_cluster(current_cluster);
        
    end
    
end

function plot_single_cluster(current_cluster)

    handles = gethand; 

    % clear the figure panel
    cla;
    colormap jet;
    legend('off');
    colorbar('delete');
    axis auto;
    
    session_data = retr('sessionData'); % all data
    gates        = retr('gates');
    gate_context = retr('gateContext'); % indices currently selected
    if isempty(gate_context)
        return;
    end
    
    selected_channels = get(handles.lstChannels, 'Value');
    channel_names     = retr('channelNames');
    selected_gates    = get(handles.lstGates, 'Value');
    gate_names        = gates(selected_gates, 1);
    
    %find cluster channel
    cluster_channel      = get(handles.lstClusterChannels, 'Value');
    if isempty(cluster_channel) || ~isDiscrete(cluster_channel) || length(unique(session_data(gate_context, cluster_channel))) > 500, %if channels is empty or channel is not discrete or the number of clusters is too large
%         [s,v] = listdlg('PromptString','Select a cluster channel (must be discrete):',...
%                 'SelectionMode','single',...
%                 'InitialValue', initTimeSelection,...
%                 'ListString',channel_names);
        if ~v
            return;
        elseif ~isDiscrete(cluster_channel)
            return;
        elseif length(unique(cluster_channel)) > 500,
            return;
        end
        cluster_channel = s;
    end
    
    %finding single cluster plot type (1:marker density, 2:heat map, 3:graph)
    plot_type      = get(handles.popSingleClusterPlotType, 'Value');
   
    if plot_type == 1,  %marker density
        
        
        %L2_dist = L2_dist_heatmap_data(session_data(gate_context, selected_channels), session_data(gate_context, cluster_channel));
        %L2_dist = L2_dist.*sign(L2_dist);   %changing all negative values to positive values (ie changing to just L2)

        SPR_distances = SPR_dist_heatmap(session_data(gate_context, selected_channels), session_data(gate_context, cluster_channel), 3);
        SPR_distances = SPR_distances.*sign(SPR_distances); %changing all negative values to positive so that markers where distributions are most different can be found
        
        %[~,sig_markers_rank] = sort(L2_dist(current_cluster,:), 'descend');
        [~,sig_markers_rank] = sort(SPR_distances(current_cluster,:), 'descend');
        
%         if sum(sig_markers) == 0,   %if no significant markers are found for this cluster
%             fprintf(strcat('No significant markers were found for cluster:',mat2str(current_cluster),'\n'));
%         end

        top = str2num(get(handles.txtSingleClusterDensTop, 'String'));   %number of markers to show, choosen by user (default = 9)
        
        if top > length(selected_channels),
            top = length(selected_channels);
        end

        %figuring out how to divide plot into subplots to fit all significant markers
        y_number = round(sqrt(top));
        if y_number < sqrt(top),    %if rounded square root is less than square root
            x_number = y_number+1;
        else
            x_number = y_number;
        end

    
        whole_fs = zeros(top,100);
        whole_ixs = zeros(top, 100);
        cluster_fs = zeros(top,100);
        cluster_ixs = zeros(top, 100);
        
        for k=1:top;   %looping through significant marker channels
            marker = sig_markers_rank(k);
            [whole_fs(k,:),whole_ixs(k,:)] = ksdensity(session_data(gate_context,selected_channels(marker)));    %density of all clusters
            [cluster_fs(k,:),cluster_ixs(k,:)] = ksdensity(session_data(session_data(gate_context,cluster_channel) == current_cluster,selected_channels(marker)));  %density of curent cluster

        end
        
        ymax = max(max(max(whole_fs)),max(max(cluster_fs)))+0.2;
            
        p = 1;  %indicator of which subplot plot should be plotted in

        for k=1:top,
            marker = sig_markers_rank(k);
            subplot(x_number, y_number,p), area(whole_ixs(k,:),whole_fs(k,:),'FaceColor',[0.65,0.65,0.65], 'linestyle', 'none');
            ylim([0,ymax]); 
            hold;
            subplot(x_number, y_number,p), plot(cluster_ixs(k,:),cluster_fs(k,:), 'r', 'linewidt', 2);
            set(gca, 'xtick', []);
            set(gca,'ytick',[]);
            
            set(title(channel_names(selected_channels(marker))));

            
            p = p+1;
        end
  
    
    elseif plot_type == 2,  %Heat map
        
        color = get(handles.popSampleHeatColor, 'Value');
        
        if color == 1,
            color_map = interpolate_colormap(redbluecmap, 64);
        elseif color == 2,
            color_map = color_map_creator('rg');
        end
        
        %calculating input to plots
        %hl_L2_dist = L2_dist_heatmap_data(session_data(gate_context, selected_channels), session_data(gate_context, cluster_channel));
        SPR_distances = SPR_dist_heatmap(session_data(gate_context, selected_channels), session_data(gate_context, cluster_channel), 3);

        heat_min = min(min(SPR_distances));
        heat_max = max(max(SPR_distances));
        
        %scaling colors so that zero will always be black
        if abs(heat_min) > heat_max,
            heat_max = abs(heat_min);
        elseif heat_min < 0,
            heat_min = -heat_max;
        end
        
%         L2_dist = hl_L2_dist.*sign(hl_L2_dist);   %changing all negative values to positive values (ie changing to just L2)
        %L2_density_acceptance = L2_marker_significance(L2_dist);    %finding markers whos L2 distance is higher than threshold for significance
        
        marker_std = std(session_data(session_data(gate_context,cluster_channel) == current_cluster, selected_channels)) ./ std(session_data(gate_context, selected_channels));
        
        
        %plotting
        subplot(40,1,1:4), imagesc(SPR_distances(current_cluster,:), [heat_min,heat_max]);
        colormap(color_map);
        set(gca, 'Ytick', []);
        set(gca, 'Xtick', []);
        %colorbar;
        %cb=colorbar;
        %set(cb,'YTick',[heat_min,heat_max]) 
        set(title('Distance between cluster distribution and whole sample distribution'));

        freezeColors;
        cbfreeze(colorbar);

%         subplot(11,1,2), imagesc(L2_dist(current_cluster,:));
%         colormap(color_map_creator('grey_scale_wb'));
%         colorbar
%         %cb=colorbar;
%         %set(cb,'YTick',[0,1])        
%         set(gca, 'Ytick', []);
%         set(gca, 'xtick', linspace(1.5,length(selected_channels)-0.5, length(selected_channels)-1));
%         set(gca, 'XTickLabel', []);
%         colorbar('hide');
        
%         freezeColors;
%         cbfreeze(colorbar);
        
        subplot(40,1,7:10), imagesc(marker_std);
        colormap(color_map_creator('grey_scale_bw'));
        colorbar;
        %cb=colorbar;
        %set(cb,'YTick',[min(marker_std),max(marker_std)])
        set(gca, 'Ytick', []);
        xticklabel_rotate([1:length(selected_channels)],90,channel_names(selected_channels));
        set(gca, 'Xtick', []);
        set(title('Ratio of within and across cluster standard deviations'));

        
    elseif plot_type == 3,  %graph
        
        %computing data for graph
        compute_dremi_graph(session_data(session_data(gate_context,cluster_channel) == current_cluster, selected_channels), current_cluster, channel_names(selected_channels));
        fprintf('Graph has been computed an will be visible in new window\n');
    end
    
end


function createMeta
    handles = gethand;
    
    gates        = retr('gates');
    sessionData  = retr('sessionData');
    gate_context = retr('gateContext');
    selected_channels = get(handles.lstChannels,'Value');
    gate_names        = get(handles.lstGates, 'String');
    selected_gates = get(handles.lstGates, 'Value');
    [~, channel_names] = getSelectedIndices(selected_gates);
    
    
    %find channels that might be cluster channels => discrete channels
    
    params = retr('metaParams');
    if (isempty(params))
        params = [];
        params.neighbors = 2;
        params.metric = 'euclidean';
        params.selected_gates = gate_names;
        params.all_channels = channel_names;
        params.cluster_channel = 1;
        %params.all_channels = numel(channel_names');
    end
    
    params = create_metaGUI('gates', gate_names, 'params', params);
    if (isempty(params))
        return;
    end
    put('metaParams', params);


    try
        
        %finding indeces of selected samples (currently assuming that no overlap occures between indeces)
        inds = {};
        for i=1:length(selected_gates),
            inds{i} = find(ismember(gate_context, gates{selected_gates(i),2}));
        end
        
        [~, ~, meta_cluster_channel] = LOUVAIN_meta_clusters(sessionData(gate_context,selected_channels), sessionData(gate_context,params.cluster_channel), inds, params.neighbors, params.metric);
        addChannels({'meta_clusters'}, meta_cluster_channel, gate_context);

    catch e
        uiwait(msgbox(sprintf('Finding meta clusters failed: %s', e.message),...
            'Error','error'));  
        return;        
    end
    
end


function plot_meta_clusters(cluster_channel)
    %fprintf('Meta clusters not implemented yet\n');
    
    handles = gethand; 

    % clear the figure panel
    cla;
    colormap jet;
    legend('off');
    colorbar('delete');
    axis auto;
        
    % show controls
    if ~isCtrlPressed
        set(handles.pnlMetaClusterControls, 'Visible', 'on');
    end

	hPlot = subplot(1,1,1,'Parent',handles.pnlPlotFigure);
    box on;

    session_data = retr('sessionData'); % all data
    gates        = retr('gates');
    gate_context = retr('gateContext'); % indices currently selected
    if isempty(gate_context)
        return;
    end
    
    selected_channels = get(handles.lstChannels, 'Value');
    channel_names     = retr('channelNames');
    selected_gates    = get(handles.lstGates, 'Value');
    gate_names        = gates(selected_gates, 1);
    
    initClusterSelection = size(session_data, 2);
    
    %find meta cluster channel
    meta_channel      = get(handles.lstMetaClusterChannels, 'Value');
    if isempty(meta_channel) || ~isDiscrete(meta_channel) || length(unique(session_data(gate_context, meta_channel))) > 500, %if channels is empty or channel is not discrete or the number of clusters is too large
%         [s,v] = listdlg('PromptString','Select a cluster channel:',...
%                 'SelectionMode','single',...
%                 'InitialValue', initClusterSelection,...
%                 'ListString',channel_names);
        if ~v
            return;
        elseif ~isDiscrete(meta_channel)
            return;
        elseif length(unique(meta_channel)) > 500,
            return;
        end
        meta_channel = s;
    end
    
    %find cluster channel
    cluster_channel      = get(handles.lstBasisOfMetaChannel, 'Value');
    if isempty(cluster_channel) || ~isDiscrete(cluster_channel) || length(unique(session_data(gate_context, cluster_channel))) > 500, %if channels is empty or channel is not discrete or the number of clusters is too large
%         [s,v] = listdlg('PromptString','Select a cluster channel:',...
%                 'SelectionMode','single',...
%                 'InitialValue', initClusterSelection,...
%                 'ListString',channel_names);
        if ~v
            return;
        elseif ~isDiscrete(cluster_channel)
            return;
        elseif length(unique(cluster_channel)) > 500,
            return;
        end
        cluster_channel = s;
    end
    
    
    %finding indeces of selected samples (currently assuming that no overlap occures between indeces)
    inds = {};
    
    for i=1:length(selected_gates),
        inds{i} = find(ismember(gate_context, gates{selected_gates(i),2}));
    end
    
    %finding cluster centroids and mapping between clusters and meta clusters
    meta_cluster_channel = session_data(gate_context, meta_channel);
    [centroids, cluster_mapping] = get_centroids(session_data(gate_context, selected_channels), inds, session_data(gate_context, cluster_channel), meta_cluster_channel); 
    
    %finding plot type
    plot_types = get(handles.popMetaPlotOptions, 'String');
    plot_type = plot_types(get(handles.popMetaPlotOptions, 'Value'));
   
    if strcmp(plot_type,'Meta cluster distribution'),
        
        cells_in_meta = zeros(max(cluster_mapping(:,1)),size(selected_gates,2));
        
        for i=1:length(selected_gates),  %looping through samples
            for j=1:max(cluster_mapping(:,1)),    %looping through meta clusters
                                
                cells_in_meta(j,i) = sum(meta_cluster_channel(inds{i}) == j)/size(inds{i},1); %percentage of cells in sample belonging to meta cluster
            end
        end

        %plotting stacked bar graph
        subplot(20, 1, 1:15), bar(cells_in_meta', 0.9, 'stack');
        ylim([0,1]);
        legend(strread(num2str(1:size(unique(meta_cluster_channel),1)),'%s'),'location','eastoutside');
        title('Percentage of cells beloning to each meta cluster');
        ylabel('% cells');
        %xticklabel_rotate([1:length(selected_gates)],90,gate_names);  
        xticklabel_rotate([1:length(selected_gates)],90,gate_names);  

        colormap(distinguishable_colors(size(unique(meta_cluster_channel),1)));
        
    elseif strcmp(plot_type,'tSNE map'),
        
        if ~isCtrlPressed
            set(handles.pnlMetaTsneColor, 'Visible', 'on');
        end
       
        %testing that no centroids are NaN (can happen if subsample does
        %not contain any cells from a specific cluster)
        unique_meta_clusters = unique(cluster_mapping(:,1));
        cluster_mapping(isnan(centroids(:,1)),:) = [];    %removing rows that correspond to clusters where centroid is NaN
        
        %testing if any meta clusters get dropped (can happen eg happen if meta cluster only contains one cluster and this contains NaN in centroids
        if (length(unique_meta_clusters) ~= length(unique(cluster_mapping(:,1)))),
            fprintf('\nA meta cluster is not plotted becaus no cells in selected gates belong to this meta cluster\n\n');
        end
        
        centroids(isnan(centroids(:,1)),:) = [];  %removing rows with NaNs in centroids
        
        dot_size = 500/max(cluster_mapping(:,5)) * cluster_mapping(:,5)+5;  %rescaling size so they make sense size wise
        
        %finding previous tSNE or calculating tSNE
        tSNE_out = retr('tsneParams');
        if (isempty(tSNE_out)),
            if (size(centroids, 1) > 10)
                tSNE_out = fast_tsne(centroids, 50, 10);    %running tSNE on centroids
            else
                tSNE_out = tsne(centroids, [], 2, size(centroids,2));  
            end
        end
        put('tsneParams', tSNE_out);
        
        %finding color channel
        color_chan = get(handles.lstTsneColorBy,'Value');        
        
        %testing if color_chan is discrete (= maybe cluster channel) or continous (maybe marker channel)
        %if isDiscrete(color_chan) && length(unique(session_data(gate_context, color_chan))) < 500,
        if color_chan == meta_channel,
            tsne_col = cluster_mapping(:,1);
            scatter_by_point(tSNE_out(:,1), tSNE_out(:,2), tsne_col, dot_size); %plotting

        else
            
            %finding marker means for marker selected
            inds = {};
    
            for i=1:length(selected_gates),
                inds{i} = find(ismember(gate_context, gates{selected_gates(i),2}));
            end
            
            tsne_col = zeros(size(centroids,1),1);
            for i=1:size(centroids,1),
                sub_data = session_data(inds{cluster_mapping(i,3)},:);
                tsne_col(i) = mean(sub_data(sub_data(:,cluster_channel) == cluster_mapping(i,4), color_chan));
                %tsne_col(cluster_mapping(:,1) == i) = median(session_data(session_data(gate_context,cluster_channel) == i,color_chan));
            end
            
            %making 0.95 quantile most red color and 0.05 quantile most blue color to compensate for outliers
            tsne_col(tsne_col < quantile(tsne_col, 0.05)) = quantile(tsne_col,0.05);
            tsne_col(tsne_col > quantile(tsne_col, 0.95)) = quantile(tsne_col, 0.95);
            
            scatter(tSNE_out(:,1),tSNE_out(:,2), dot_size, tsne_col, 'fill');
            colorbar;

        end
        
        xlabel('tSNE1');
        ylabel('tSNE2');
        
        %tSNE_out = fast_tsne(centroids, 50, 10);
        %tSNE_plot(tSNE_out, 'metaclusters',meta_cluster_list, out_dir, dot_size);
        %scatter_by_point(tSNE_out(:,1), tSNE_out(:,2), color_variable, dot_size)

              
        
    elseif strcmp(plot_type,'Single meta cluster'),
        
        % show controls
        if ~isCtrlPressed
            set(handles.pnlSingleMetaControls, 'Visible', 'on');
        end
        
        unique_meta = unique(meta_cluster_channel);
        set(handles.lstSingleMetaCluster, 'String', unique_meta);
        
        %find selected meta cluster
        current_meta = unique_meta(get(handles.lstSingleMetaCluster, 'Value'));

        %generating heatmap data for all clusters for all time points
        heatmap_data = [];
        ylabs = [];
        
        for i=1:length(selected_gates),   %looping through gates
            heat = SPR_dist_heatmap(session_data(gate_context(inds{i}), selected_channels), session_data(gate_context(inds{i}), cluster_channel), 3);
            heatmap_data = vertcat(heatmap_data,heat(cluster_mapping(cluster_mapping(:,1) == current_meta & cluster_mapping(:,3) == i,4),:));
            ylabs = vertcat(ylabs,strcat(gate_names(i),'(',strread(num2str(cluster_mapping(cluster_mapping(:,1) == current_meta & cluster_mapping(:,3) == i,4)'), '%s'),')'));
        end
        
        %finding min and max values for heat map
        heat_min = min(min(heatmap_data));
        heat_max = max(max(heatmap_data));
        
        %scaling colors so that zero will always be black
        if abs(heat_min) > heat_max,
            heat_max = abs(heat_min);
        elseif heat_min < 0,
            heat_min = -heat_max;
        end
        
        %defining color map
%         color = get(handles.popSampleHeatColor, 'Value');
%         
%         if color == 1,
%             color_map = interpolate_colormap(redbluecmap, 64);
%         elseif color == 2,
%             color_map = color_map_creator('rg');
%         end
        color_map = interpolate_colormap(redbluecmap, 64);
        
        %plot heat map
        subplot(10,1,1:9), imagesc(heatmap_data, [heat_min,heat_max]);
        %set(title(strcat('Sample ', mat2str(p), ': L2 + higher/lower than median')));
        set(gca, 'Position', [0.18,0.27,0.8,0.66]);
        set(gca, 'ytick', 1:length(ylabs));
        set(gca, 'Yticklabel', ylabs);
        set(gca, 'xtick', []);
        colormap(color_map);
        colorbar;
        xticklabel_rotate([1:length(selected_channels)],90,channel_names(selected_channels));
        
        
        %fprintf('Function not implemented yet\n');
    end
    
end

function generateConditionGates

    %function generates new gates
    
    handles = gethand;
    
    gates        = retr('gates');
    sessionData  = retr('sessionData');
    gate_context = retr('gateContext');
    selected_channels = get(handles.lstChannels,'Value');
    gate_names        = get(handles.lstGates, 'String');
    selected_gates = get(handles.lstGates, 'Value');
    channel_names = retr('channelNames');
    
    %test that selected channel is discrete => can be cluster channel
    %takes too long for large datasets
%     if ~isDiscrete(sessionData(gate_context,selected_channels)) || length(unique(session_data(gate_context, selected_channels))) > 500
%         fprintf('Selected channel is not discrete, condition gates could not be generated\n')
%         return;
%     end
    
    
    %find binary channels
    binary_channels = zeros(length(channel_names),1);
    
    for i=1:length(channel_names)
%         if sum(sessionData(gate_context, i) == 0 | sessionData(gate_context, i) == 1) == size(sessionData(gate_context,:),1)
%             binary_channels(i) = 1;
%         end
        if all(sessionData(gate_context, i) == 0 | sessionData(gate_context, i) == 1)
            binary_channels(i) = 1;
        end
    end
    
    if sum(binary_channels) == 0
        fprintf('No binary condition channels found. Can not generate new gates\n');
        return;
    end

    %run GUI interface
    params = retr('conditionParams');
    if (isempty(params))
        params = [];
        params.permutations = 1000;
        params.condition_channels = channel_names(logical(binary_channels));
        params.inc_markers = channel_names(~logical(binary_channels));
        params.reference = channel_names(logical(binary_channels));
    end    
    put('conditionParams', params);

    params = create_ConditionGatesGUI('params', params);
    
     if (isempty(params))
        return;
     end

    %ask for path to put session data files once new gates have been computed
    pathname = uigetdir('dialog_title', 'Save new gates in folder');

    if isequal(pathname,0)   %eg if user pressed cancel
        return;
    end

    
    %generate new gates
    try
        
        
        %find indexes of channels
        binary_idx = find(binary_channels == 1);
        ref_channels = binary_idx(params.reference);
        cond_channels = binary_idx(params.condition_channels);
        cond_channels(cond_channels == ref_channels) = [];
        
        nonbinary_idx = find(binary_channels == 0);
        inc_markers = nonbinary_idx(params.inc_markers);
  
        
        %finding indecies for selected gates in sessionData(gate_context,:)
        inds = {};
        for i=1:length(selected_gates),
            inds{i} = find(ismember(gate_context, gates{selected_gates(i),2}));
        end
        
        %finding new channel names
        new_channel_names = {'cluster_ID'; 'cell_number'; 'cell_percentage'};
        new_channel_names = vertcat(new_channel_names, strcat('basal_',channel_names(inc_markers))');
        
        selected_data = sessionData(gate_context,:);   %variable only containing data form selected gates
        
        for i=1:length(inds)    %looping through samples
            tic;

            gates = retr('gates');  %gates is changed for each loop, should be updated
            
            
            data = selected_data(inds{i},:);    %finding data belonging to sample
            unique_clusters = unique(sessionData(inds{i}, selected_channels));
            
            new_gate = zeros(length(unique_clusters), 3+length(inc_markers)*(length(cond_channels)+1)); %one new gate for each sample
            
            for j=1:length(unique_clusters) %looping through clusters in sample
            
                %adding cluster ID
                new_gate(j,1) = unique_clusters(j);
                
                %adding number of cells in sample belonging to cluster
                new_gate(j,2) = size(data(data(:,selected_channels) == unique_clusters(j),:),1);
                
                %adding percentage of cells in sample belonging to cluster
                new_gate(j,3) = size(data(data(:,selected_channels) == unique_clusters(j),:),1) / size(data,1);
                
                %tic;
                %adding mean marker levels for reference for cluster
                new_gate(j,4:3+length(inc_markers)) = mean(data(data(:,selected_channels) == unique_clusters(j) * sum(data(:,ref_channels) == 1, 2) == 1,inc_markers));
                %fprintf('Computing basal means sample %s, cluster %s. time:%g\n',gate_names{selected_gates(i)}, mat2str(unique_clusters(j)), toc);

                
                %find 
                for k=1:length(cond_channels)   %looping through each condition
                    
                    if j == 1 && i == 1 %if this is first sample and first cluster (channel names the same for all samples and all clusters)
                        new_channel_names = vertcat(new_channel_names, strcat(channel_names(cond_channels(k)),'_',channel_names(inc_markers))');
                    end
                    
                    %finding distances between condition distributions for cluster and basal
                    distances = zeros(1,length(inc_markers));
                    
                    for p=1:length(inc_markers) %looping through markers
                        cluster_distribution = data(data(:,selected_channels) == unique_clusters(j) * data(:,cond_channels(k)) == 1,inc_markers(p));
                        reference_distribution = data(data(:,selected_channels) == unique_clusters(j) * sum(data(:,ref_channels) == 1, 2) >= 1,inc_markers(p));
                        
                        %calculating distance
                        perm = str2double(params.permutations);
                        %tic;
                        [score,~,~,~,~] = SAER(reference_distribution, cluster_distribution, perm);
                        %fprintf('Computing distance condition %s, marker %s. time:%g\n',channel_names{cond_channels(k)}, channel_names{inc_markers(p)}, toc);

                        distances(p) = score;

                        
                    end
                    
                    new_gate(j, 4+length(inc_markers)*k:3+length(inc_markers)*(k+1)) = distances;
                    
                end

            end
            
            %fprintf('Computing gate for sample %s. time:%g\n',gate_names{selected_gates(i)}, toc);
            
            %adding new gate to GUI
            %tic;

%             tic;
%             cluster_gate_inds = size(sessionData, 1)+1:...
%                                 size(sessionData, 1)+size(new_gate,1);
%             fprintf('Finding new gate indecies %s.time:%g\n\n',gate_names{selected_gates(i)}, toc);
%             
%             tic;
%             sessionData(end+1:end+size(new_gate,1), :) = 0;
%             put('sessionData', sessionData);
%             fprintf('Adding to session data %s. time:%g\n\n',gate_names{selected_gates(i)}, toc);
%             
%             tic;
%             add_gate(...
%                 strcat(gate_names{selected_gates(i)},'_conditions'),...
%                 cluster_gate_inds,...
%                 {});
%             fprintf('Adding gate %s. time:%g\n\n',gate_names{selected_gates(i)}, toc);
%             
%             tic;
%             addChannels(reshape(new_channel_names, 1, numel(new_channel_names)),...
%                 new_gate,...
%                 cluster_gate_inds,...
%                 size(gates, 1)+1);
%             fprintf('Adding channels %s. time:%g\n\n',gate_names{selected_gates(i)}, toc);
            
            file_name = strcat(pathname, '/', gate_names{selected_gates(i)},'_condition.csv');
            csvwrite_with_headers(file_name,new_gate,new_channel_names);
            
            fprintf('Computing and saving gate for sample %s. time:%g\n\n',gate_names{selected_gates(i)}, toc);
            
        end 
            
    catch e
        uiwait(msgbox(sprintf('Generating condition gates failed: %s', e.message),...
            'Error','error'));  
        return;        
    end
    
end



function hPlot = plotDensity
    handles = gethand; 

    % clear the figure panel
    cla;
    colormap jet;
    legend('off');
    colorbar('delete');
    set(gca, 'CLimMode', 'auto');
    axis auto;
    
    % show controls
    if ~isCtrlPressed
        set(handles.pnlDensityControls, 'Visible', 'on');
    end

    session_data = retr('sessionData'); % all data
    gate_context = retr('gateContext'); % indices currently selected
    if isempty(gate_context)
        return;
    end
    
	hPlot = subplot(1,1,1,'Parent',handles.pnlPlotFigure);
    box on;

    selected_channels = get(handles.lstChannels, 'Value');
    channel_names     = retr('channelNames');

    set(handles.pupDensXAxis, 'Value', selected_channels(1));
    set(handles.pupDensYAxis, 'Value', selected_channels(2));

    vx = session_data(gate_context, selected_channels(1));
    vy = session_data(gate_context, selected_channels(2));

    my_plot_dens(vx, vy,0);
           
	xlabel(channel_names{selected_channels(1)}, 'Interpreter', 'none')
    ylabel(channel_names{selected_channels(2)}, 'Interpreter', 'none')
end

function lstHistGroup_Callback
    handles = gethand;
    selectedGates       = get(handles.lstHistGroupA, 'Value');
    put('selectedGates', selectedGates);
   
    plot_histograms(1);        
end

function plotEnrichment
    handles      = gethand;

    selected_channels = get(handles.lstChannels, 'Value');
    selected_gates    = get(handles.lstGates, 'Value');
    
    sessionData  = retr('sessionData');
    gates        = retr('gates');
    gateContext  = retr('gateContext');
    channelNames = retr('channelNames');
    channelNames = channelNames(selected_channels);
    
    %-- clear the figure panel
    if ~isCtrlPressed
        hidePlotControls;
    end

    cla;
    colormap jet;
    legend('off');
    colorbar('delete');
    set(gca, 'CLimMode', 'auto');
    axis auto;
    title('');
    
    if (isempty(gateContext) || isempty(selected_channels)) return; end
    
    hPlot = subplot(1,1,1,'Parent',handles.pnlPlotFigure);

    %-- end clear

    
    nplots = numel(selected_gates);
    nrows = floor(nplots/6)+1;
    ncols = ceil(nplots/nrows);
    
    for i=1:nplots
        
        subplot(nrows, ncols, i);
        
            colormap(genColorMap('rw', 20));
            
            data = sessionData(gates{selected_gates(i), 2}, selected_channels);
            means = mean(data)';
            [sorted_means, idx] = sort(means, 1, 'descend');
            imagesc(sorted_means, [0 max(sorted_means)]);
            set(gca,'ytick',1:length(channelNames))
            set(gca,'yticklabel', channelNames(idx), 'fontsize', 9)
            colorbar;
            title(gates{selected_gates(i), 1});
            
    end
end

% There are many types of scatter plots and unfortunately each requires
% a different matlab function:
% 2D or 3D X Colored by Gate or by channel or Diff of two channels
function hPlot=plotScatter 

    handles      = gethand;
    sessionData  = retr('sessionData');
    gates        = retr('gates');
    gateContext  = retr('gateContext');
    channelNames = retr('channelNames');

    selected_channels = retr('selectedChannels');
    nSelChannels     = numel(selected_channels);
    
    selected_gates = get(handles.lstGates, 'Value');
    % if no channel selected => select all channels
    if (numel(selected_gates)==0 || nSelChannels==0) 
        selected_gates = 1:size(gates, 1);
    end

    % make sure axis popups are visible
    if ~isCtrlPressed
        hidePlotControls;
        set(handles.pnlScatterControls, 'Visible','on');
    end
    
    % clear the figure panel
    cla;
    colormap jet;
    legend('off');
    colorbar('delete');
    set(gca, 'CLimMode', 'auto');
    set(gca,'XTickLabelMode','auto');
    
    axis auto;
    title('');
    
    if (isempty(gateContext) || isempty(selected_channels)) return; end
    
    hPlot = subplot(1,1,1,'Parent',handles.pnlPlotFigure);
    box on;
    
    nCH1 = selected_channels(1);
    nCH2 = selected_channels(2);
    nCH3 = 0;
    set(handles.pupXAxis, 'Value', nCH1);
    set(handles.pupYAxis, 'Value', nCH2);
    if (nSelChannels == 3) 
        nCH3 = selected_channels(3);
        set(handles.pupZAxis, 'Value', nCH3+1);
    else
        set(handles.pupZAxis, 'Value', 1);        
    end
    
%     plot_knn(handles, nSelChannels, sessionData, gateContext, nCH1, nCH2, nCH3);
%     hold on;
   
    nChColors = get(handles.pupColorBy, 'Value')-1;
    nplots = numel(nChColors);
    nrows = floor(nplots/3)+1;
    ncols = ceil(nplots/nrows);
    
    for i=1:nplots
        nChColor = nChColors(i);
        subplot(nrows, ncols, i);
        
    % color by gate
    if (nChColor == 0)
        
        % add filter
        intIndices = getIntIndices;

        matColors = distinguishable_colors(max(selected_gates));
        
        hold on;
        
        aggregateinds = [];
        vColor = [];
        for gi=selected_gates
            currGate = gates{gi, 2};
            if (~isempty(intIndices))
                currGate = intersect(intIndices, currGate);
            end
            vColor = [vColor; ones(numel(currGate), 1)*gi];
            aggregateinds = [aggregateinds; currGate(:)];
        end
        vX = sessionData(aggregateinds, nCH1);
        vY = sessionData(aggregateinds, nCH2);

            if (nSelChannels == 3)
                vZ = sessionData(aggregateinds, nCH3);
            else
                if get(handles.chkRandomLayers, 'Value')
                    vZ = rand(1, numel(vX))*10;
                else
                    vZ = vColor;
                end
            end

        % plot 
        if get(handles.chkGradientCommunities, 'Value')
            clr = jet(numel(selected_gates));
        else
            clr = distinguishable_colors(numel(selected_gates));
%             clr(2, :) = []; % remove red
        end
        
            vColor_discrete = vColor;
            colors = unique(vColor)';
            for ci=1:numel(colors);
                vColor_discrete(vColor==colors(ci)) = ci;
            end

        
%         if (~get(handles.chkDensity, 'Value'))
            myplotclr(vX, vY, vZ, vColor_discrete, 'o', clr, [min(vColor_discrete), max(vColor_discrete)], nSelChannels == 3);
%         end
        box on;
        colorbar off;

%         vSize = ones(numel(gateContext), 1);
%         vSize(gates{2, 2}) = 2;
%         vSize = mynormalize(vSize, 99.99)*1000;
%         myscatter_by_point(vX, vY, vZ, vColor, vSize, true);
%         box on;
%         colorbar off;

        if get(handles.chkLegend, 'Value') && numel(selected_gates) > 1
            l = legend(remove_repeating_strings(gates(selected_gates, 1)), 'Interpreter', 'none');
            if (numel(selected_gates) > 6)
                set(l, 'Location','NorthEastOutside');
            end
        end
        fixLegendMarkerSize;
        hold off;
    % color by channel
    else
        % scatter all selected gates and color by a channel
        if ((~retr('diff') || numel(selected_gates) ~= 2) || nSelChannels ==3)
            vX = sessionData(gateContext, nCH1);
            vY = sessionData(gateContext, nCH2);
            vColor = sessionData(gateContext, nChColor);

            unqValues = unique(vColor);
            if (numel(unqValues) > 100 || ~isDiscrete(nChColor) )
 
                if get(handles.chkOutlier, 'Value')
                    clim = adjustClimsToUserScaling([prctile(vColor, 0.3)  prctile(vColor, 99.7)]);
                else
                    clim = adjustClimsToUserScaling([min(vColor) max(vColor)]);
                end

                if (nSelChannels == 3)
                    vZ = sessionData(gateContext, nCH3);
                else
                    vZ = rand(1, numel(vX))*(clim(2)-clim(1))+clim(1);
                end

                
%         if (~get(handles.chkDensity, 'Value'))
                myplotclr(vX, vY, vZ, vColor, 'o', colormap, clim, nSelChannels == 3);   
%         end
                colorbar;
                caxis(clim);
                
                
%                   vSize = sessionData(gateContext, end);
%                  vSize = mynormalize(vSize, 99.99)*1000;
%                   myscatter_by_point(vX, vY, vColor, vSize, false);

            else % color by community.
                
                if get(handles.chkGradientCommunities, 'Value')
                    clr = jet(numel(unqValues));
                else
                    clr = distinguishable_colors(numel(unqValues));
                end
                
                    
                
                if (nSelChannels == 3)
                    vZ = sessionData(gateContext, nCH3);
                else
                    if get(handles.chkRandomLayers, 'Value')
                        vZ = rand(1, numel(vX))*10;
                    else
                        vZ = vColor;
                    end
                end

            vColor_discrete = vColor;
            colors = unique(vColor)';
            for ci=1:numel(colors);
                vColor_discrete(vColor==colors(ci)) = ci;
            end
            
            % plot 
%             if (~get(handles.chkDensity, 'Value'))
            myplotclr(vX, vY, vZ, vColor_discrete, 'o', clr, [min(vColor_discrete), max(vColor_discrete)], false)
%             end
            colorbar off;
            hl=legend(cellfun(@(n)(num2str(n)), num2cell(unique(vColor)), 'UniformOutput', false));                    

%               gscatter(vX, vY ,vColor, clr, '.', get(0, 'DefaultLineMarkerSize'), doleg);

%                  vSize = sessionData(gateContext, end);
%                  vSize = mynormalize(vSize, 99.99)*1000;
%                   myscatter_by_point(vX, vY, vColor, vSize, true);

            fixLegendMarkerSize;
            end
            
        % color the difference of a channel between two gates
        else  
            basal_gate = gates{selected_gates(1), 2};
            drug_gate  = gates{selected_gates(2), 2};
            
            % filter indices if user selected to intersect gates
            if get(handles.btnIntersect, 'Value')
                [inds, ~] = getSelectedIndices(get(handles.lstIntGates, 'Value'));
                if (~isempty(inds))
                    basal_gate = intersect(inds, basal_gate);
                    drug_gate  = intersect(inds, basal_gate);
                end
            end
            
            vX = sessionData(basal_gate, nCH1);
            vY = sessionData(basal_gate, nCH2);
            vZ = sessionData(basal_gate, nChColor);

            vA = sessionData(drug_gate, nCH1);
            vB = sessionData(drug_gate, nCH2);
            vC = sessionData(drug_gate, nChColor);


            if get(handles.rdbDiffBox, 'Value')
                thresh_val = get(handles.sldrThreshold, 'Value');
                box_size_val = get(handles.sldrBox, 'Value');
                thresh = round(thresh_val);
                nbins = round(30-(box_size_val*2));

                selected_ind = union(basal_gate, drug_gate);
                XandA = sessionData(selected_ind, nCH1);
                YandB = sessionData(selected_ind, nCH2);

                %# bin centers (integers)
                xbins = linspace(floor(min(XandA)),ceil(max(XandA)),nbins);
                ybins = linspace(floor(min(YandB)),ceil(max(YandB)),nbins);
                xNumBins = numel(xbins); yNumBins = numel(ybins);

                %# map X/Y values to bin indices
                Xi = round( interp1(xbins, 1:xNumBins, vX, 'linear', 'extrap') );
                Yi = round( interp1(ybins, 1:yNumBins, vY, 'linear', 'extrap') );

                %# count number of elements in each bin
                numPtsXY = accumarray([Yi(:) Xi(:)], 1, [yNumBins xNumBins]);
                sumXY = accumarray([Yi(:) Xi(:)], vZ, [yNumBins xNumBins]);

                % TODO take care of thresh
                XY = sumXY./numPtsXY;
    %             XY(find(numPtsXY<thresh)) = 0;

                %# map A/B values to bin indices
                Ai = round( interp1(xbins, 1:xNumBins, vA, 'linear', 'extrap') );
                Bi = round( interp1(ybins, 1:yNumBins, vB, 'linear', 'extrap') );

                %# count number of elements in each bin
                numPtsAB = accumarray([Bi(:) Ai(:)], 1, [yNumBins xNumBins]);
                sumAB = accumarray([Bi(:) Ai(:)], vC, [yNumBins xNumBins]);

                % TODO take care of thresh
                AB = sumAB./numPtsAB;

                diff = AB-XY;
                diff(find(numPtsAB<thresh)) = NaN;
                diff(find(numPtsXY<thresh)) = NaN;
                    
                %# plot 2D histogram
                nanimagesc(diff, genColorMap('rwb', 20), 'zeroc', 1, 'xbins', xbins, 'ybins', ybins);
                colorbar; set(gca,'YDir','normal'); hold off;
%                 hold on, plot(vX, vY, 'b.', 'MarkerSize', 10)
%                 plot(vA, vB, 'g.', 'MarkerSize', 4), 
%                 hold off;
            elseif get(handles.rdbDiffVoronoi, 'Value')

                thresh = round(get(handles.sldrDiffVoronoi, 'Value'));
                tic
                [v,c]=voronoin([vX vY]);
                disp(sprintf('Voronoi cells computed: %gs',toc));

                tic
                hold on;
                vColors = NaN(length(c), 1);
                for ci = 1:length(c) 
                    if all(c{ci}~=1)   % If at least one of the indices is 1, 
                                      % then it is an open region and we can't 
                                      % patch that.
                        in = find(inpoly([vA vB], v(c{ci},:)));
                        if (length(in)>=thresh) 
                            vColors(ci) = mean(vC(in))-vZ(ci);
                            patch(v(c{ci},1),v(c{ci},2),vColors(ci)); % use color i.
                        end
                    end
                end
                
                [cmap clim] = setColorbar(vColors);

                disp(sprintf('Voronoi cells colored: %gs',toc));
                colorbar;
                set(gca,'CLim',clim);
                hold on, plot(vX, vY, 'b.', 'MarkerSize',4)
                plot(vA, vB, 'g.', 'MarkerSize',4);
                hold off;
            else
                plot_KNN_diff(nChColor);
                return;
            end
        end
    end
    
    
    xlabel(channelNames{nCH1}, 'Interpreter', 'none');
    ylabel(channelNames{nCH2}, 'Interpreter', 'none');
    
    if nChColor > 0
        title(channelNames{nChColor},'Interpreter', 'none');
    elseif nplots > 1
        title('Color by gates','Interpreter', 'none');
    end
    
    if (nSelChannels == 3)
        zlabel(channelNames{nCH3});
        view(3);
    end

    end
    % --- add density plot for 2D
    if (nSelChannels == 2)
        gating = 'on';
        if nplots > 1
            gating = 'off';
        end
        enableGating(handles, gating);
        view(2);
        h = zlabel(' ');
        set(h, 'String', '');

        if (get(handles.chkDensity, 'Value'))
            [~, density, x, y] = kde2d(sessionData(gateContext, [nCH1 nCH2]), 512);
            hold on;
%             cmap = jet;
%             cmap(1, :) = [1, 1, 1];
%             colormap(cmap);
            contour(x, y, density, 512);
            hold off;
        end
        if (get(handles.chkLog, 'Value'))
            set(gca, 'XScale', 'log');
            set(gca, 'YScale', 'log');
        else
            set(gca, 'XScale', 'linear');
            set(gca, 'YScale', 'linear');
        end  
        
        if (get(handles.chkFixAxis, 'Value'))
            try
                clim_xmin = str2num(get(handles.txtMinXAxis, 'String'));
                clim_xmax = str2num(get(handles.txtMaxXAxis, 'String'));
                clim_ymin = str2num(get(handles.txtMinYAxis, 'String'));
                clim_ymax = str2num(get(handles.txtMaxYAxis, 'String'));
                clims = [clim_xmin clim_xmax clim_ymin clim_ymax];
                axis(clims);
            catch
                setStatus('Warning: please enter numerical values for manual colorbar scaling');
            end    
        end
   end
    tic
    % --- add KNN plot
    plot_knn(handles, nSelChannels, sessionData, gateContext, nCH1, nCH2, nCH3);
end

function fixLegendMarkerSize

    handles = gethand;
    
    % fix legend marker sizes
    if get(handles.chkLegend, 'Value')
        l = findobj(gcf,'tag','legend');
        a=get(l,'children');

        % a(i%3) corresponds to the marker object
        % i starts from 5 from Matlab 2013 and on.
        tic
        try 
            for k=1:size(a, 1)
                if strcmpi(get(a(k),'Type'), 'line') 
                    set(a(k),'markersize',6);
                end
            end
        catch e
            disp 'error while trying to increase legend marker size';
        end
        toc;
    end
end

function enableGating(handles, state)
    set(handles.btnGatePoly, 'Enable', state);
    set(handles.btnGateRect, 'Enable', state);   
end

function plot_KNN_diff(color_ch)
    handles = gethand;
    session_data = retr('sessionData');
    gates = retr('gates');
    
    selected_gates = get(handles.lstGates, 'Value');
    channel_names = gates{selected_gates(1), 3};
    [basal_gate ~] = getSelectedIndices(selected_gates(1));
    [stim_gate ~] = getSelectedIndices(selected_gates(2));
    
    if get(handles.btnIntersect, 'Value')
        intInds = getIntIndices;
        basal_gate = intersect(basal_gate, intInds);
        stim_gate = intersect(stim_gate, intInds);
    end
    
    % get knn space channels
    knn_space = get(handles.lstKNNSpace, 'Value');
    
    % get color by channel
%     color_ch = get(handles.pupColorBy, 'Value') - 1;
    
    % get thresh and K
    K = round(get(handles.sldrDiffKNN, 'Value'));
    try 
        thresh = str2double(get(handles.txtDiffKNNThresh, 'Value'));
    catch
        setStatus(['KNN distance threshold must be a number, current value ' thresh]);
        thresh = 0;
    end
    
    % get knn space data from both and presentation dims
    b_data = session_data(basal_gate, knn_space);
    s_data = session_data(stim_gate, knn_space);
    b_exp = session_data(basal_gate, color_ch);
    s_exp = session_data(stim_gate, color_ch);
    
    presentation_dim = get(handles.lstKNNPresentationSpace, 'Value');
    
    if numel(b_exp) > 25000 || numel(s_exp) > 25000
        selection = questdlg(sprintf('You''ve chosen to run KNN on %g points and this will take some time, are you sure you want to go ahead?', numel(b_exp)+numel(s_exp)) ,...
             'KNN on many points',...
             'Yes','No','No');
         if strcmp(selection,'No')
            return
         end
    end
    
    % if rdbknnsourceboth
    if get(handles.rdbKNNSourceBoth, 'Value')
        [b_colors s_colors] = FCmap(b_data, s_data, b_exp, s_exp, 'thresh', thresh, 'K', K, 'stim', 1);
        expression = [b_colors ; s_colors];
        presentation_space = session_data([basal_gate(:); stim_gate(:)], presentation_dim);
    else % only basal
        [b_colors ~] = FCmap(b_data, s_data, b_exp, s_exp, 'thresh', thresh, 'K', K);
        expression = b_colors;
        presentation_space = session_data(basal_gate, presentation_dim);
    end
    [cmap clim] = setColorbar(expression);
    scatter(presentation_space(:, 1), presentation_space(:, 2), [],expression, '.');

    colorbar;
    caxis(clim);
    set(gca, 'Color', [.7 .7 .7]);
    
    xlabel(channel_names{presentation_dim(1)}, 'Interpreter', 'none');
    ylabel(channel_names{presentation_dim(2)}, 'Interpreter', 'none');
    
    gate_names = (gates(selected_gates, 1));
    gate_names = remove_repeating_strings(gate_names);
    str_title = sprintf('%s minus %s :: %s', gate_names{2}, gate_names{1}, channel_names{color_ch});
    title(str_title,'Interpreter', 'none');

    enableGating(handles, 'on');
    view(2);
    h = zlabel(' ');
    set(h, 'String', '');
end

function clims=adjustClimsToUserScaling(clims)
    handles = gethand;
    if get(handles.chkClim, 'Value')
        try
            clim_min = str2num(get(handles.txtMinColor, 'String'));
            clim_max = str2num(get(handles.txtMaxColor, 'String'));
            clims = [clim_min clim_max];
        catch
            setStatus('Warning: please enter numerical values for manual colorbar scaling');
        end
    end
end

function [cmap clim]=setColorbar(vColors)
	nancolor = [0.7 0.7 0.7];
    set(gca, 'Color', nancolor);
    cmap = genColorMap('rwb', 20);
    nclevel = size(cmap,1); %nclevel is an odd number
    cmap(2:end+1,:) = cmap;    
    cmap(1,:) = nancolor; %set NaN to gray, NaN is consider min
    m = nanmin(nanmin(vColors));
    M = nanmax(nanmax(vColors));        

    if isnan(m)
        m = 0;
    end
    if isnan(M)
        M = 1;
    end
    if M > 0 && m < 0
        midcolor = (nclevel+1)/2 + 1; %index of mid-level color in colormap
        halfcolor = (nclevel-1)/2;
        if M > -m
            i = ceil(-m/M*halfcolor);
            cmap = [cmap(1,:); cmap(midcolor-i:end,:)];
        else
            i = ceil(-M/m*halfcolor);
            cmap = [cmap(1,:); cmap(2:midcolor+i,:)];
        end
        nclevel = size(cmap,1) - 1;

%         display(sprintf('min: %g', m));
%         display(sprintf('colorbar min: %g', m-(M-m)/(nclevel-1)));
%         display(sprintf('max: %g', M));
%         display(sprintf('num colors: %g', size(cmap,1)));
        step_size = (M-m)/(nclevel-1);
        [~, white_idx] = ismember([1 1 1], cmap, 'rows'); % get white idx
       
        clim_min = -(white_idx-0.5)*step_size;
        clim_max = (size(cmap,1)- white_idx+0.5)*step_size;
%         display(sprintf('cmin: %g', clim_min));
%         display(sprintf('cmax: %g', clim_max));
%         display(sprintf('-------'));
    else
        step_size = (M-m)/(nclevel-1);
        clim_min = m-step_size;
        clim_max = M;
    end
    colormap(cmap);
    clim = [clim_min clim_max];
end

function plot_knn(handles, nSelChannels, sessionData, gateContext, nCH1, nCH2, nCH3)
	% draw KNN graph
	k_values = get(handles.pupKNN, 'String');
    selection    = get(handles.pupKNN, 'Value');
    k_neighbors = str2num(k_values{selection});
    if (k_neighbors)
        if (nSelChannels == 2)
            mat = sessionData(gateContext, [nCH1 nCH2]);
        else 
            mat = sessionData(gateContext, [nCH1 nCH2 nCH3]);
        end
        IDX = knnsearch(mat,mat, 'K', k_neighbors+1);
        for i=1:size(IDX,1)
            for j=2:size(IDX,2)
                pts = [mat(i, :);  mat(IDX(i, j), :)];
                if (nSelChannels == 2)
                    line(pts(:, 1), pts(:,2));
                else
                    line(pts(:, 1), pts(:,2), pts(:,3));
                end
            end
        end
    end

end

function k_sparse = create_k_sparse(mat, k_neighbors)
    % compute adjacency mat
    adj = pdist2(mat,mat);
	adj = abs(adj);
    
    % create a sparse KNN matrix
	IDX = knnsearch(adj,adj, 'K', k_neighbors);
    j_sparse_mat = IDX';
    i_col = (1:size(adj,1))';
    i_sparse_mat = i_col*ones(1, k_neighbors);
    i_sparse_mat = i_sparse_mat';
    k_sparse = sparse(i_sparse_mat(:), j_sparse_mat(:), 1);


%     IDX = knnsearch(mat,mat, 'K', k_neighbors);
%     k_sparse = sparse(size(mat,1),size(mat,1));
%     k_sparse(IDX) = 1;
end

function runTSNE(normalize)
    ndims = 2; % fast tsne is only implemented for 2 dims.
    handles = gethand;
    
    selected_channels = get(handles.lstChannels,'Value');
            
    sessionData  = retr('sessionData');
    gate_context = retr('gateContext');

    MAX_TSNE = 1000000;
    
    if (numel(gate_context) > MAX_TSNE)
        setStatus(sprintf('Cannot run tSNE locally on more than %g points. Please subsample first.', MAX_TSNE));
        return;
    end
    
    
    hwaitbar = waitbar(0,'BH-SNE started. Refer to Command window or std-o for progress...');
% 	waitbar(0,hwaitbar, sprintf('Begin tSNE on %g points.\n Please be patient and check standard output for progress updates', numel(gate_context)));
%     map = compute_mapping(sessionData(gate_context, selected_channels), 'tSNE', ndims);
    tic;
    data = sessionData(gate_context, selected_channels);
        
    selection = questdlg('Compute over normalized data?' ,...
                     'Normalize',...
                     'Yes','No','Yes');
    if strcmp(selection,'Yes')
        data = mynormalize(data, 99);
    end

%     new_selected_channels = selected_channels(prctile(data, 99) >= 2.2);
%     set(handles.lstChannels,'Value', new_selected_channels);
%     data(:, prctile(data, 99) < 1.2) = [];
      
    map = fast_tsne(data, 110);
%     map = tsne(data, [], ndims);
    disp(sprintf('map generated in %g m', toc/60));
    setStatus('Done tSNE calc. Setting results to gates');

% 	waitbar(2/n, hwaitbar);
    new_channel_names = cell(1, ndims);
    for i=1:numel(new_channel_names)
        new_channel_names{i} = sprintf('bh-SNE%g', i);
    end

    addChannels(new_channel_names, map, gate_context);

	waitbar(1,hwaitbar, 'Done.');
    
    close(hwaitbar);
	setStatus('Done.');
    
end

function runWanderlust
    handles = gethand;
    
    gates        = retr('gates');
    sessionData  = retr('sessionData');
    gate_context = retr('gateContext');
    selected_channels = get(handles.lstChannels,'Value');
    gate_names        = get(handles.lstGates, 'String');
    
    data = sessionData(gate_context, selected_channels);
    
    params = retr('wanderlustParams');
    if (isempty(params) || ~isfield(params, 'l'))
        params = [];
        params.l = 15;
        params.k = 8;
        params.num_landmarks = 100;
        params.num_graphs = 10;
        params.metric = 'euclidean';
        params.normalize = true;
        params.selected_gate = numel(gate_names);
    else
        params.selected_gate = min(params.selected_gate,numel(gate_names));
    end
    
    params = wanderlustGUI('gates', gate_names, 'params', params);
    if (isempty(params) || ~isfield(params, 'l'))
        return;
    end
    put('wanderlustParams', params);
	s = find(ismember(gate_context,gates{params.selected_gate, 2}));
    
    if (isempty(s))
        uiwait(msgbox('Selected starting points are not in the gate context (selected gates)',...
            'Invalid starting gate','error'));               
        return;
    end
    si = randsample(1:numel(s), 1);
    params.s = s(si);
    
    try
        % normalize data is asked
        if (params.normalize)
            data = data-repmat(prctile(data, 1, 1), size(data,1),1);
            data = data./repmat(prctile((data), 99, 1),size(data,1),1);

            data(data > 1) = 1;
            data(data < 0) = 0;
            data(isinf(data)) = 0;
        end
              
        % cycler params
        params.band_sample = true;
        params.voting_scheme = 'exponential';
        params.flock_landmarks = 2;
%         params.plot_landmark_paths = true;
%         params.plot_data = data(:,1:2);
        
        % compute trajectory
        G = wanderlust(data,params);
                
        % Save result
        addChannels({'cct'}, mean(G.T, 1)', gate_context);
        
    catch e
        fprintf('Cycler failed: %s', getReport(e,'extended'));
        uiwait(msgbox(sprintf('Cycler failed: %s', getReport(e,'basic')),...
            'Error','error'));  
        return;        
    end
    
end

function runLouvain(handles)
    if askuser('Would you like to phenograph each gate separately?')
    	phenoEach();
    else 
        handles = gethand;

        session_data  = retr('sessionData');
        gate_context = retr('gateContext');
        selected_channels = get(handles.lstChannels,'Value');

        % ask user for cluster count to recognize
        k_neigh = inputdlg('Enter number of neigh: ',...
                              'cluster each', 1, {num2str(floor(sqrt(numel(gate_context)/2)))}); 

        if (isempty(k_neigh) || str2num(k_neigh{1}) == 0)
            return;
        end

        k_neigh = str2num(k_neigh{1});

        data = session_data(gate_context, selected_channels);

        [IDX, ~] = phenograph(data, k_neigh);
        addChannels({sprintf('pheno%g', k_neigh)}, IDX);
    end
      % -- Lovain implementation. 
%     handles = gethand;
%     
%     sessionData  = retr('sessionData');
%     gate_context = retr('gateContext');
%     selected_channels = get(handles.lstChannels,'Value');
%     
%     louvain_res = louvain(sessionData(gate_context,:), selected_channels,...
%         'channelnames', retr('channelNames'));
%     
%     if (isempty(louvain_res)) 
%         % user cancelled
%         return;
%     end
%     
%     % add new louvain channels
%     mods = louvain_res{2};
% 
%     new_channel_names = cell(1, min (numel(mods), numel(louvain_res{1})));
%     for i=1:numel(new_channel_names)
%         new_channel_names{i} = sprintf('louvain_%g', mods(i));
%     end
%     
%     addChannels(new_channel_names, cell2mat(louvain_res{1}));
%     
	setStatus('Done.');
end

% new_channels_names - 1XN Cell array of strings
% new_data           - MXN data matrix
% opt_gate_context   - opt: SessionData indices. Default is gateContext
% opt_gates          - opt: the gate index to add the channel name on to.
function addChannels(new_channel_names, new_data, opt_gate_context, opt_gates)
    handles = gethand;

    % set variables
    sessionData  = retr('sessionData');
    gates        = retr('gates');    

    if (exist('opt_gate_context','var'))
        gate_context = opt_gate_context;
    else
        gate_context = retr('gateContext');    
    end
    
    if (exist('opt_gates','var'))
        selected_gates = opt_gates;
    else
        selected_gates = get(handles.lstGates, 'Value');
        if isempty(selected_gates)
            selected_gates = 1:size(gates, 1);
        end
        
        % filter indices if user selected to intersect gates
        if get(handles.btnIntersect, 'Value')
            selected_int_gates = get(handles.lstIntGates, 'Value');
            [gate_indices channel_names] = getSelectedIndices(selected_gates);
            [gate_int_indices channel_int_names] = getSelectedIndices(selected_int_gates);
            % check if one group is contained in the other
            if isempty(setdiff(gate_int_indices, gate_indices))
                selected_gates = selected_int_gates;
            else
                msgbox('Your are using intersect mode so SightOf does not know which gates to add the resulting channels to. By default, when the intersecting group is not contained in the main selected gates group, the channels are added to all the main selected gates. ','Channels added to selected gates though content is only added to the intersection.','warn');
            end
        end    

    end
    
    % add necessary channels to the selected gates
    defined_channels = cellfun(@(x)numel(x), gates(selected_gates, 3), 'uniformoutput', true);
    undef_channel_ind = max(defined_channels)+1;
    
    if (size(sessionData,2)-undef_channel_ind >= 0) && ...
        any(~any(sessionData(gate_context, undef_channel_ind:end)))
        
        % find a streak the same width of new_data of empty columns
        d = diff([false any(sessionData(gate_context, undef_channel_ind:end)) == 0 ones(1, size(new_data, 2)) false]);
        p = find(d==1);
        m = find(d==-1);
        lr = find(m-p>=size(new_data, 2));
        last_def_channel = undef_channel_ind - 1 + (p(lr) - 1);
    else
        last_def_channel = size(sessionData,2);
    end
        
    for i=selected_gates
        
        % add new tsne channel names to gate
        channel_names = gates{i, 3};
        if (last_def_channel-numel(channel_names) > 0)
            % add blank\placeholder channel names
            for j=numel(channel_names)+1:last_def_channel
                channel_names{j} = 'cyt_placeholder_tmp';
            end
        end
        channel_names(end+1:end+numel(new_channel_names)) = new_channel_names;
        gates{i, 3} = channel_names;
    end
    
    n_new_columns = size(new_data, 2) - (size(sessionData,2) - last_def_channel);
    
    % extend session data
    if (n_new_columns > 0)
        new_columns = zeros(size(sessionData, 1), n_new_columns);
        sessionData = [sessionData new_columns];
    end
    
    % set new data to session
    sessionData(gate_context, last_def_channel+1:last_def_channel+size(new_data, 2)) = new_data;    

    
    put('sessionData', sessionData);
    put('gates', gates);
    
    lstGates_Callback;
end

function addGate(name, inds)
    gates         = retr('gates'); % all gates (names\indices) in cell array
    channel_names = retr('channelNames');
    
    gates{end+1, 1} = name;
    gates{end, 2}   = inds;
    gates{end, 3}   = channel_names;

    put('gates', gates);
end

function cluster(isKmeans)
    handles = gethand;
    session_data = retr('sessionData');
    gate_context = retr('gateContext');
    selected_channels = get(handles.lstChannels, 'Value');

    % ask user for cluster count to recognize
    n_clusters = inputdlg('Enter number of clusters: ',...
                          'K-Means', 1, {num2str(floor(sqrt(numel(gate_context)/2)))}); 
    
    if (isempty(n_clusters) || str2num(n_clusters{1}) == 0)
        return;
    end

    k_clusters = str2num(n_clusters{1});
    data = session_data(gate_context, selected_channels);
        
    if (isKmeans)        
        % run K-Means
        clust_alg = 'KMeans';
    	IDX = kmeans(data, k_clusters,'Display', 'iter', 'EmptyAction', 'singleton');
    else 
        % run EMGM
        clust_alg = 'EMGM';
        [IDX, model, llh] = emgm(session_data(gate_context, selected_channels)', k_clusters);
        IDX = IDX';
    end
    
    % add results to GUI
    addChannels({sprintf('%s%g',clust_alg, k_clusters)}, IDX);
end

function runAffinityPropegation
    handles = gethand;
    session_data = retr('sessionData');
    gate_context = retr('gateContext');
    selected_channels = get(handles.lstChannels, 'Value');
    tic;
%     N = numel(gate_context);
%     x=session_data(gate_context, selected_channels);
%     M=N*N-N; s=zeros(M,3); 
%     j=1;
%     for i=1:N
%       for k=[1:i-1,i+1:N]
%         s(j,1)=i; s(j,2)=k; s(j,3)=-sum((x(i,:)-x(k,:)).^2);
%         j=j+1;
%       end;
%     end;
%     p=median(s(:,3)); % Set preference to median similarity

    A = session_data(gate_context, selected_channels)';
    A = bsxfun(@minus,A,mean(A,2));
    S = full(dot(A,A,1));
    fprintf('dot product: %g', toc); 
    tic
    D = bsxfun(@plus,S',S)-full(2*(A'*A));
    fprintf('created pairwise dist matrix: %g', toc); 
    p = median(D(:));
    s = D;
    
    % reformat imput mat to a N^2 by 3 values
    N = size(D, 1);
    [X Y] = meshgrid(1:N, 1:N);
    s = [X(:), Y(:) ,D(:)];
    k = 1:N;
    s(N*(k-1)+k, :) = []; % clear diag
    
    [IDX,netsim,dpsim,expref]=apclusterSparse(s,p,'plot', 'maxits', 40);
    fprintf('Number of clusters: %d\n',length(unique(IDX)));
    fprintf('Fitness (net similarity): %f\n',netsim);

    % add results to GUI
    addChannels({'affinity_propegation'}, IDX);

    x = D;
    figure; % Make a figures showing the data and the clusters
    for i=unique(IDX)'
      ii=find(IDX==i); h=plot(x(ii,1),x(ii,2),'o'); hold on;
      col=rand(1,3); set(h,'Color',col,'MarkerFaceColor',col);
      xi1=x(i,1)*ones(size(ii)); xi2=x(i,2)*ones(size(ii)); 
      line([x(ii,1),xi1]',[x(ii,2),xi2]','Color',col);
    end;
    axis equal tight;
end
% --- Executes on button press in chkDensity.
function chkDensity_Callback(hObject, eventdata, handles)
    plotScatter;
end

% --- Executes on selection change in pupPlotType.
function pupPlotType_Callback(~, ~, ~)
end

% --- Executes on selection change in pupScatterPlotType.
function pupScatterPlotType_Callback(~, ~, ~)
end

% --------------------------------------------------------------------
function btnGate_ClickedCallback(~, ~, ~, isPoly)
    handles = gethand;
    sessionData  = retr('sessionData');
    gateContext   = retr('gateContext');
    selectedChannels = retr('selectedChannels');
    selGates = get(handles.lstGates, 'Value');
    
    CH1 = selectedChannels(1);
    CH2 = selectedChannels(2);

    setStatus('Waiting: Double click on node when finished');
    gate = tsne_gate(gcf, sessionData(gateContext, [CH1 CH2]), 1, isPoly);
    setStatus(sprintf('you have gated: %g data points', numel(gate))); 

    if size(gate, 2) > 0
        % create new gate
        created=createNewGate(gateContext(gate), retr('channelNames'));

        % if more than one gate was selected - add a channel 'gate source'
        gates = retr('gates');        
        if created && numel(selGates) > 1
            sessionData     = retr('sessionData');
            nDataLength     = size(sessionData,1);
            indGateSource   = size(gates, 1);
            gateSources     = zeros(nDataLength, 1);

            % create a new channel for gate course index.
            for selGate=selGates
                gateSources(gates{selGate, 2}) = selGate;
            end
            
            addChannels({'gate_source'}, gateSources, 1:nDataLength, indGateSource);
        end
    end
end

function btnGateVal_Callback
    handles = gethand;
    sessionData  = retr('sessionData');
    gateContext   = retr('gateContext');
    selectedChannels = retr('selectedChannels');
    CH1 = selectedChannels(1);
    CH2 = selectedChannels(2);

    setStatus('Waiting: Double click on node when finished');
    gate = tsne_gate(gcf, sessionData(gateContext, [CH1 CH2]), 1, isPoly);
    setStatus(sprintf('you have gated: %g data points', numel(gate))); 

    if size(gate, 2) > 0         
        createNewGate(gateContext(gate), retr('channelNames'));

        gates         = retr('gates');
%         set(handles.lstGates, 'String', gates(:, 1));
%         set(handles.lstIntGates, 'String', gates(:, 1));
    end
    
end
% --- Executes on selection change in Scatter's axis combo boxes.
function pupAxis_Callback(~, ~, ~)
    handles = gethand;

    if (get(handles.pupZAxis, 'Value')==1) 
        set(handles.lstChannels, 'Value', [get(handles.pupXAxis, 'Value') get(handles.pupYAxis, 'Value')]);
        lstChannels_Callback;
        put('selectedChannels', [get(handles.pupXAxis, 'Value') get(handles.pupYAxis, 'Value')]);
    else
        set(handles.lstChannels, 'Value', [get(handles.pupXAxis, 'Value') get(handles.pupYAxis, 'Value') get(handles.pupZAxis, 'Value')-1]);
        lstChannels_Callback;
        put('selectedChannels', [get(handles.pupXAxis, 'Value') get(handles.pupYAxis, 'Value') (get(handles.pupZAxis, 'Value')-1)]);
    end
    plotScatter;
end

% --- Executes on selection change in pupXAxis.
function pupDensAxis_Callback
    handles = gethand;
    
    selectedAxis = [get(handles.pupDensXAxis, 'Value') get(handles.pupDensYAxis, 'Value')];
    
    set(handles.lstChannels, 'Value', selectedAxis);
    put('selectedChannels', selectedAxis);

    plotDensity;
end

% --- Executes on selection change in pupMarkerSize.
function pupMarkerSize_Callback(hObject, ~, ~)
    handles = gethand;
    marker_values = get(hObject, 'String');
    selection    = get(hObject, 'Value');
    marker_size = str2num(marker_values{selection});
    set(0,'DefaultLineMarkerSize',marker_size);
    put('markerSize', marker_size);
    plotChannels_Callback(0, 0, handles);
end

% --- Executes on button press in btnIntersect.
function btnIntersect_Callback(hObject, ~, ~)
    
    if get(hObject,'Value')
    	% Toggle button is pressed-take appropriate action
        set(hObject, 'cdata', brighten(double(get(hObject, 'cdata')), double(-0.5)));
        
        handles = gethand;
        
        pos = get(handles.lstGates, 'Position');
        pos(4) = pos(4)/2;
        pos(2) = pos(2) + pos(4);
        set(handles.lstGates, 'Position', pos);
  
        set(handles.lstIntGates, 'Visible', 'on');
    else
    	% Toggle button is not pressed-take appropriate action
        set(hObject, 'cdata', brighten(double(get(hObject, 'cdata')), double(double(0.5))));

        handles = gethand;
        pos = get(handles.lstGates, 'Position');
        pos(2) = pos(2) - pos(4);
        pos(4) = pos(4)*2;
        set(handles.lstGates, 'Position', pos);
  
        set(handles.lstIntGates, 'Visible', 'off');
    end
    
    lstGates_Callback;
end

function refreshGates
    handles = gethand;
    gates = retr('gates');
    n = size(gates, 1);
    
    %reset selected index if necessary
    resetListSelection(handles.lstGates, n);
    resetListSelection(handles.lstIntGates, n);
    resetListSelection(handles.lstHistGroupA, n);
    resetListSelection(handles.lstHistGroupB, n);
    
    % replace gate names
    gate_names = gates(:, 1);
    set(handles.lstGates, 'String', gate_names);
    set(handles.lstIntGates, 'String', gate_names);
    set(handles.lstHistGroupA, 'String', gate_names);
    set(handles.lstHistGroupB, 'String', gate_names);
end

function resetListSelection(lstPanel, maxSelection)
    selectedGate = get(lstPanel, 'Value');
    if (selectedGate > maxSelection)
        set(lstPanel, 'Value', 1);
    end
end
function toggleRightHist()
    handles = gethand;
    toggleHist(handles.btnHistIntRight, 'right');
end

function toggleLeftHist()
    handles = gethand;
    toggleHist(handles.btnHistIntLeft, 'left');
end

%toggles image on toggle hObject button
function toggleHist(hObject, dir)  
    if get(hObject,'Value')
       placeIcon(hObject,  ['arrow' dir '.png']);
    else
       placeIcon(hObject,  ['arrow' dir '_inac.png']);
    end
    plot_histograms(1);
end

% --- Executes on button press in btnSaveGate.
function btnSaveGate_Callback(~, ~, ~)
    handles = gethand;
    
    channelNames = retr('channelNames');
    sessionData = retr('sessionData');
    selGateInd   = retr('gateContext');
    gates = retr('gates');
    selectedGates = get(handles.lstGates, 'value');
    
    if isempty(gates)
        return;
    end
    
    % if only a single gate is selected and it has a filename associated
    % already then we'll suggest to overwrite original file name.
    if ((numel(selectedGates) == 1 ) && ~isempty(gates{selectedGates, 4}))
        [filename,pathname,~] = uiputfile('*.fcs','Save Gate', gates{selectedGates, 4});
    else
        [filename,pathname,~] = uiputfile('*.fcs','Save Gate');
    end

    if isequal(filename,0) || isequal(pathname,0)
        return;
    end
    
    % remove placeholders
    columns = find(~ismember(channelNames,'cyt_placeholder_tmp')==1);
    
    % write the gate\s to the specified FCS filename.
    fca_writefcs([pathname filename], sessionData(selGateInd, columns), channelNames(columns),channelNames(columns));
    
    % associate selected gate to filename incase a single gate was selected
    if (numel(selectedGates) == 1)
        gates{selectedGates, 4} = [pathname filename];
        put('gates', gates);
    end
end

% --- Executes on button press in btnLoadGates.
function btnLoadGates_Callback(~, ~, ~)
    OpenMenuItem_Callback
end

% --- Executes on button press in btnSplitCellCycle.
function btnSplitCellCycle(~,~,~)

    session_data = retr('sessionData');
    gate_context = retr('gateContext');
    featurenames = retr('channelNames');  
    
    
%     gates = retr('gates');
%     handles = gethand;
%     selected_gates  = get(handles.lstGates, 'Value'); % currently selected    
%     ris = [];
%     for i=selected_gates
% 
%         gate_context = gates{i, 2};
    
    % evaluate phase using data and specific channels (selected by name)
    data = session_data(gate_context, :);
    try 
        [G1, S, G2, M0, M1] = cyclerclassification_v2(featurenames, data);
    catch e
         uiwait(msgbox(sprintf('Automatic classification and cleaning failed: %s', getReport(e,'basic')),...
        'Error','error'));  
        disp(getReport(e,'extended'));     
    end
 
    % add indicator channels
    delG1 = zeros(size(gate_context)); delG1(G1) = 1;
    delS  = zeros(size(gate_context)); delS(S)  = 1;
    delG2 = zeros(size(gate_context)); delG2(G2) = 1;
    delM0 = zeros(size(gate_context)); delM0(M0) = 1;
    delM1 = zeros(size(gate_context)); delM1(M1) = 1;
    addChannels({'G1','S','G2','M0','M1'},...
                [delG1(:), delS(:), delG2(:), delM0(:), delM1(:)]);

%         ris(end+1) = ri;

    % split cells to gates by phase
    addGate('G1', gate_context(G1));
    addGate('S', gate_context(S));
    addGate('G2', gate_context(G2));
    addGate('M0', gate_context(M0));
    addGate('M1', gate_context(M1));
    interphase = union(union(gate_context(G1),gate_context(S)), gate_context(G2));
    addGate('Interphase', interphase);
    addGate('Start', randsample(gate_context(G1), 1));
    
%     end

%     mean(ris(:))
%     ris(:)
end

% --- Executes on button press in btnRemoveGate.
function btnReorderGate_Callback(dir)
    handles = gethand;

    gates = retr('gates');
    selectedGates = get(handles.lstGates, 'Value');
    
    if isempty(selectedGates) || isempty(gates) ||... 
        (dir < 0 && min(selectedGates) == 1)    || ...
        (dir > 0 && max(selectedGates) == size(gates, 1))
        return;
    end

    % move gates in direction given
    order = 1:size(gates, 1);
    prev   = setdiff(selectedGates+dir, selectedGates);
    newpos = setdiff(selectedGates, selectedGates+dir);    
    order([newpos selectedGates+dir]) = order([prev selectedGates]);
    gates = gates(order, :);
    
    %save order
    put('gates', gates);
    
    % update gui
%     set(handles.lstGates, 'String', gates(:, 1));
%     set(handles.lstIntGates, 'String', gates(:, 1));
    if (~isempty(gates))
        set(handles.lstGates, 'Value', selectedGates+dir);
    end
    
    lstGates_Callback;

end

% --- Executes on button press in btnRemoveGate.
function btnRemoveGate_Callback(~, ~, ~)
    handles = gethand; 

    % remove gates from gate list
    gates = retr('gates');
    remGates = get(handles.lstGates, 'Value');
    gates(remGates, :) = [];
    put('gates', gates);
%     set(handles.lstGates, 'String', gates(:, 1));
%     set(handles.lstIntGates, 'String', gates(:, 1)); 
    % we keep all of the session data because if we remove it we'll have
    % to update each gate's indices as the indices will have changed. 
%     % cleanup session data values that have no pointers from gates
% 	sessionData = retr('sessionData');
%     keepIndices = [];
%     for i=1:size(gates, 1)
%         keepIndices = union(gates{i, 2}, keepIndices);
%     end
%     sessionData = sessionData(keepIndices, :);
%     put('sessionData', sessionData)

    % update gui
    if (~isempty(gates))
        set(handles.lstGates, 'Value', max(remGates(1)-1, 1));
    else
        set(handles.lstGates, 'Value', 0);
        put('sessionData', []);
    end
    
    lstGates_Callback;
end

% --------------------------------------------------------------------
function showGettingStarted
    open Cycler_WorkflowManual.pdf;
end

function cmiAbout_CreateFcn(hObject, ~, ~)

end

function cmiTransformGate_Callback(~, ~, ~)
    handles = gethand;
    sessionData = retr('sessionData');
    gates = retr('gates');
    selGates = get(handles.lstGates, 'Value');
    
    for i=selGates
        if isempty(gates{i,4}) % Gate is associated with a file.
            uiwait(msgbox('Transformation\cleanup is only done for gates that are associated with a specific file.','Cannot perform transformation on selected files.','modal'));
            return;
        end
    end
    
    % assume all gates have the same channels.
    for i=selGates
        if (~exist('cofactor', 'var'))
            [cofactor selectedChannels isDnagate isSaveout isOriginal isPrefix strPrefix] = ...
                Preprocess('channelNames', gates{i, 3});
            if cofactor == 0
                return;
            end
            hwaitbar = waitbar(0, 'preprocess ...')
        end

        waitbar((find(selGates == i)*2-1)/(numel(selGates)*2), hwaitbar, sprintf('transforming %s', gates{i,1}));
        setStatus(sprintf('transforming %s', gates{i,1}));

        gateData = sessionData(gates{i, 2}, :);

        % Transform
        if (~isempty(selectedChannels))
            gateData(:, selectedChannels) = asinh( gateData(:, selectedChannels) ./ cofactor );
        end
        
        % save gate data locally (this will also affect other gates
        % pointing to same data. should not be a problem since we're
        % only preprocessing file associated data
        sessionData(gates{i, 2}, :) = gateData;

        % DNA gate
        if (isDnagate)
            debris_threshold = 0.9

            chDNA = cellstrfnd(gates{i, 3}, 'DNA');
            DNA = gateData(:, chDNA);

            dna_gm = gmdistribution.fit(DNA, 2);
            [~,j] = max(sum(dna_gm.mu,2));
            P = dna_gm.posterior(DNA);
            removeDNA = find(P(:,j) < debris_threshold);
%                 gateData(removeDNA,:) = []; TODO find a way to delete the
%                 data. This line would mess up the indices for events
%                 after the deletion.
            gates{i, 2}(removeDNA) = []; % <-- remove the indices from the gate
        end

        put('sessionData', sessionData);
        put('gates', gates);

        if (isSaveout)
            savePath = gates{i,4};
            if isPrefix
                [PATH,NAME,EXT] = fileparts(savePath);
                NAME = [strPrefix NAME];
                savePath = fullfile(PATH, [NAME EXT]);
            end

            waitbar(find(selGates == i)/numel(selGates), hwaitbar, sprintf('Saving %s to:\n %s', gates{i,1}, savePath));
            setStatus(sprintf('saving %s to %s', gates{i,1}, savePath));
            chNames = gates{i, 3};
            fca_writefcs(savePath,...
                gateData(:, 1:numel(chNames)),...
                chNames,chNames);
            if isPrefix
                gates{i,4} = savePath;
                gates{i,1} = NAME;
%                 set(handles.lstGates, 'String', gates(:, 1));
%                 set(handles.lstIntGates, 'String', gates(:, 1));
                put('gates', gates);
            end               
        end
    end  
    
    if (~isempty(selGates)) 
        setStatus(sprintf('Done.'));
        close(hwaitbar);
    end
end

% --------------------------------------------------------------------
function cmiRenameGate_Callback(~, ~, ~)
    handles = gethand;
    gates = retr('gates');
    selected_gates = get(handles.lstGates, 'Value');
    
    for i=selected_gates
        new_name = inputdlg(sprintf('Enter a new gate name for %s: ',...
            gates{i, 1}), 'Gate Name', 1, gates(i, 1)); 
        if (~isempty(new_name))
            gates(i, 1) = new_name;
        end
    end
    
    put('gates', gates);
%     set(handles.lstGates, 'String', gates(:, 1));
%     set(handles.lstIntGates, 'String', gates(:, 1));

end

% --------------------------------------------------------------------
function cmiRenameRegexp_Callback(hObject, eventdata, handles)
    handles = gethand;
    regexps = retr('regexps');
    
    selected_regexps = get(handles.lstRegexps, 'Value');
    if (isempty(selected_regexps))
        return;
    end
    
    reg_index = selected_regexps(1);
    new_name = inputdlg(sprintf('Enter new name for %s: ',...
        regexps{reg_index, 1}), 'Reg Exp Name', 1, regexps(reg_index, 1)); 
    if (~isempty(new_name))
        regexps(reg_index, 1) = new_name;
    end
    
    put('regexps', regexps);
    set(handles.lstRegexps, 'String', regexps(:, 1));
end

function setStatus(sStatus)
    try 
        handles = gethand;
        set(handles.lblStatus, 'String', sStatus);
    catch err
        disp err;
    end
end

function tsne_each_gate(ndims)
handles = gethand;
selected_gates = get(handles.lstGates, 'Value');
for i=selected_gates
    set(handles.lstGates, 'Value', i);
    lstGates_Callback;
    runTSNE(ndims);
end
end

function locateNonExpressivePoints(thresh)
    handles = gethand;
    selected_channels = get(handles.lstChannels, 'Value');
    selected_gates = get(handles.lstGates, 'Value');
    gateContext = retr('gateContext');
    sessionData = retr('sessionData');
    
    
    v = zeros(1, numel(gateContext));
    tic;
    vMax = max(sessionData(gateContext, selected_channels), [], 2);
    v(vMax < thresh) = 1;
    toc
    addChannels({sprintf('below_%g', thresh)}, v(:));
end

function runTsnePerGate(nDims)
	handles = gethand;
    selected_gates = get(handles.lstGates, 'Value');
    
    for i=selected_gates
        set(handles.lstGates, 'Value', i);
        lstGates_Callback;
        runTSNE(nDims);
    end

end

% --------------------------------------------------------------------
function cmiRenameChannel_Callback(hObject, eventdata, handles)
    handles = gethand;
    gates = retr('gates');
	channel_names = retr('channelNames');
    selected_channels = get(handles.lstChannels, 'Value');
    selected_gates = get(handles.lstGates, 'Value');
    
    if isempty(selected_gates)
        selected_gates = 1:size(gates, 1);
    end
    
    for i=selected_channels
        new_name = inputdlg(...
                    sprintf('Enter a new channel name for %s: ', channel_names{i}),...
                    'Channel Rename', 1, channel_names(i)); 
        if (~isempty(new_name))
            % rename channel is all selected gates
            channel_names(i) = new_name;
            for j=selected_gates
                gate_channels = gates{j, 3};
                gate_channels{i} = new_name{:};
                gates{j, 3} = gate_channels;
            end
        end
    end
    
    put('gates', gates);
%     set(handles.lstGates, 'String', gates(:, 1));
%     set(handles.lstIntGates, 'String', gates(:, 1));

    lstGates_Callback;
end

% --------------------------------------------------------------------
function cmnGates_Callback(hObject, eventdata, handles)
end

% --------------------------------------------------------------------
function cmnChannels_Callback(hObject, eventdata, handles)
    try
    handles = gethand;
    state = 'off';

    selected_channel = get(handles.lstChannels, 'Value');
    if (numel(selected_channel)==1) %&& isDiscrete(selected_channel))
        state = 'on';
        if isDiscrete(selected_channel)
            set(handles.cmiSplitToGates, 'Label', 'Split to Gates');
        else
            set(handles.cmiSplitToGates, 'Label', 'Slice to Gates...');
        end
    end
    
    set(handles.cmiSplitToGates, 'Visible', state);
    set(handles.cmiGateChannel, 'Visible', state);
	set(handles.mniCommunityLabels, 'Visible', state);
    % right click exception
    catch e
    end
end



function cmiSplitToGates_Callback(~, ~, ~)
	handles = gethand;

    selected_channel = get(handles.lstChannels, 'Value');
    selected_gates    = get(handles.lstGates, 'Value');
    gateContext  = retr('gateContext');
    sessionData  = retr('sessionData');
    gates        = retr('gates');
    channel_names = retr('channelNames');
    if isDiscrete(selected_channel)
        for selected_gate=selected_gates
            gate_data = sessionData(intersect(gates{selected_gate, 2}, gateContext), selected_channel);
            for val=unique(gate_data)'
                % --- create gate ---

                % Set Name
                gates{end+1, 1} = sprintf('%s-%s-%g',...
                    gates{selected_gate, 1},channel_names{selected_channel},val(1));

                % Set data indices (
                indices = find(sessionData(gateContext, selected_channel)==val);
                gates{end, 2} = intersect(gates{selected_gate, 2},gateContext(indices));

                % Set channel names
                gates{end, 3} = channel_names;
            end
        end
    else
        slice = slicechannels;
        if isempty(slice) || slice.num_slices < 2
            return;
        end
        
        % Compute weighted averages at each location
        x = sessionData(gateContext, selected_channel);
        locs = linspace(min(x), max(x), slice.num_slices);

        unit = locs(2)-locs(1);
        
        for selected_gate=selected_gates
            for val=locs(:)'
                % --- create gate ---

                % Set Name
                gates{end+1, 1} = sprintf('%s-%s-%g/%g',...
                    gates{selected_gate, 1},channel_names{selected_channel},find(val(1)==locs), slice.num_slices);

                % Set data indices (
                indices = find(sessionData(gateContext, selected_channel)>=(val(1)-unit*slice.overlap) &...
                    sessionData(gateContext, selected_channel)<=(val(1)+unit*slice.overlap));
                gates{end, 2} = intersect(gates{selected_gate, 2},gateContext(indices));

                % Set channel names
                gates{end, 3} = channel_names;
            end
        end
    end
    
    put('gates', gates);
% 	set(handles.lstGates, 'String', gates(:, 1));
% 	set(handles.lstIntGates, 'String', gates(:, 1));
end

function cmiGateChannel
	handles = gethand;

    selected_channel = get(handles.lstChannels, 'Value');
    selected_gates    = get(handles.lstGates, 'Value');
    gateContext  = retr('gateContext');
    gates = retr('gates');
    sessionData  = retr('sessionData');
    channel_names = retr('channelNames');
    
    x = sessionData(gateContext, selected_channel);
    params.bottom = min(x);
    params.top = max(x);
    
    gate = gatebyvalue('params', params, 'title', 'Gate Bounds');
    if isempty(gate)
        return;
    end

   % --- create gate ---

    % Set Name
    gates{end+1, 1} = sprintf('%s %g-%g',...
        channel_names{selected_channel},gate.bottom, gate.top);

    % Set data indices (
    indices = find(sessionData(gateContext, selected_channel)>=gate.bottom &...
        sessionData(gateContext, selected_channel)<=gate.top);
    gates{end, 2} = gateContext(indices);

    % Set channel names
    gates{end, 3} = channel_names;
    
	put('gates', gates);
end

%---------------------------------------------------------------------
function cmiAddChannel_Callback(~, ~, ~)
    % get filename
    % request session file
    [filename, pathname, ~] = uigetfile('*.mat', 'Load SightOff Session');
    if isequal(filename,0) || isequal(pathname,0)
        return;
    end
    
    % load file into structure
    file_contents = load(filename);
    fNames = fieldnames(file_contents);
    
    if (numel(fNames) == 1) % file has to have a single variable
        new_channel = file_contents.(fNames{1});
        gate_context = retr('gateContext');
        
        if (numel(gate_context) ~= numel(new_channel)) 
            % channel must be an array the length of gate_context
            % report error
            uiwait(msgbox(['New data channel is of size (' num2str(numel(new_channel))...
                        '). Current gate context is (' num2str(numel(gate_context)) ')'],...
                        'Invalid Channel Length','error'));               
            return;
        end
        
        % add new channel
        addChannels({fNames{1}}, new_channel(:), gate_context);

    else       
        % report error
        name_list = cellfun(@(s)([s ', ']), fNames, 'UniformOutput', false);
        name_list{end} = name_list{end}(1:end-2);
        uiwait(msgbox(['File contains more than one variable ('...
                        name_list{:} ')'],'Invalid File','error'));               
        return;    
    end
end

% --------------------------------------------------------------------
function showAbout
    aboutcyt;
end

function cmiSubsample_Callback(isSubsampleEach)
    handles = gethand;
    gateContext    = retr('gateContext');
    selected_gates = get(handles.lstGates, 'Value');
    
    

        try
            if isSubsampleEach
                txt = sprintf('You have currently selected %g points.\n\rEnter sample size: ', numel(gateContext));
            else
                txt = sprintf('You have currently selected %g gates, a total of %g points.\n\rEnter a sample size for each file: ', numel(selected_gates), numel(gateContext));
            end
            
            str_sample_size = inputdlg(...
                        sprintf(txt, numel(gateContext)),...
                        'Sample', 1, {num2str(min(2000, floor(numel(gateContext)*.20)))}); 
            if (isempty(str_sample_size) || str2num(str_sample_size{1}) == 0)
                return;
            end
            sample_size = str2num(str_sample_size{1});
            
        % while inserting non-legal number, an error msg will pop up
        catch e
            errordlg('Please anter a valid number','Wrong Input');
            return; 
        end
    
    if isSubsampleEach
        prefix = inputdlg('Enter a prefix for sample names',...
                'Gate name prefix', 1, {'sample_'});
            
        if isempty(prefix)
            return;
        end
        
        gates = retr('gates');
        channel_names = retr('channelNames');
        
        for i=selected_gates
            rand_sample = randsample(intersect(gates{i, 2}, gateContext), sample_size);
            gate_channel_names = gates{i, 3};
            if numel(channel_names) > numel(gate_channel_names)
                gate_channel_names = channel_names;
            end
            createNewGate(rand_sample, gate_channel_names, {sprintf('%s%s', prefix{1}, gates{i, 1})});
            gates = retr('gates');
%             set(handles.lstGates, 'String', gates(:, 1));
%             set(handles.lstIntGates, 'String', gates(:, 1));
        end
    else 
        rand_sample = randsample(gateContext, sample_size);
        createNewGate(rand_sample, retr('channelNames'));

        gates = retr('gates');
%         set(handles.lstGates, 'String', gates(:, 1));
%         set(handles.lstIntGates, 'String', gates(:, 1));

        % if more than one gate was selected - add a channel 'gate source'
        if numel(selected_gates) > 1
            sessionData = retr('sessionData');
            v = zeros(size(sessionData,1), 1);

            for j=selected_gates
                v(gates{j, 2}) = j;
            end
            put('sessionData', sessionData);
            addChannels({'gate_source'}, v(:), 1:numel(v), size(gates, 1));
        end
    end
end

% ----
% gate_indices - the indices in session_data that the gate points to
% channel_names - cell array of strings
% opt_gate_name - optional otherwise will request name from user.
function created=createNewGate(gate_indices, channel_names, opt_gate_name)
    created = false;
    gates = retr('gates');
    
    if ~exist('opt_gate_name','var')
        opt_gate_name = inputdlg('Enter new gate name: ', 'Gate Name');
    end
    
    if (~isempty(opt_gate_name)) 
        gates(end+1, 1) = opt_gate_name;
        gates{end, 2}   = gate_indices;
        gates{end, 3}   = channel_names;
        put('gates', gates);
        created = true;
    end
end

function isDiscrete=isDiscrete(nChannel)
    gateContext = retr('gateContext');
    sessionData = retr('sessionData');
    unique_values = unique(sessionData(gateContext, nChannel));
    isDiscrete = numel(gateContext)>0 && numel(unique_values) < 28;
    if ~isDiscrete
        if (unique_values == round(unique_values));
            isDiscrete = 1;
        end
    end
end

% --- Executes on button press in btnDiff.
function btnDiff_Callback(hObject, eventdata, handles)
    button_state = get(hObject,'Value');
    if button_state == get(hObject,'Max')
    	% Toggle button is pressed-take appropriate action
        cdata = get(hObject, 'cdata');
        put('diff', 1);
        set(hObject, 'cdata', brighten(cdata, -0.5));
        pos = get(handles.pnlPlotFigure, 'Position');
        pos(1) = pos(1)+0.083;
        set(handles.pnlPlotFigure, 'Position', pos);

        set(handles.pnlDiffControls, 'Visible', 'on');

    elseif button_state == get(hObject,'Min')
    	% Toggle button is not pressed-take appropriate action
        cdata = get(hObject, 'cdata');
        put('diff', 0);
        set(hObject, 'cdata', brighten(cdata, 0.5));
        pos = get(handles.pnlPlotFigure, 'Position');
        pos(1) = pos(1)-.083;
        set(handles.pnlPlotFigure, 'Position', pos);

        set(handles.pnlDiffControls, 'Visible', 'off');
    end 
    plotScatter;
end

% --- Executes during object creation, after setting all properties.
function lstGates_CreateFcn(hObject, ~, ~)
    if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
        set(hObject,'BackgroundColor','white');
    end    
end

% --- Executes during object creation, after setting all properties.
function pupXAxis_CreateFcn(hObject, ~, ~)
    if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
        set(hObject,'BackgroundColor','white');
    end
end

% --- Executes during object creation, after setting all properties.
function pupYAxis_CreateFcn(hObject, ~, ~)
    if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
        set(hObject,'BackgroundColor','white');
    end
end

% --- Executes during object creation, after setting all properties.
function pupColorBy_CreateFcn(hObject, ~, ~)
    if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
        set(hObject,'BackgroundColor','white');
    end
end

% --- Executes during object creation, after setting all properties.
function pupScatterPlotType_CreateFcn(hObject, ~, ~)
    if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
        set(hObject,'BackgroundColor','white');
    end
end

% --- Executes during object creation, after setting all properties.
function lstChannels_CreateFcn(hObject, ~, ~)
    if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
        set(hObject,'BackgroundColor','white');
    end
end

% --------------------------------------------------------------------
function mniCommunityLabels_Callback(hObject, eventdata, handles)
    channelNames = retr('channelNames');
    sessionData = retr('sessionData');
    selGateInd = retr('gateContext');
    data = sessionData(selGateInd, :);
    
    [channel2Labels] = label_communities(data, channelNames, channel2Labels);
end

function btnSwitchDensityAxis_Callback
    handles = gethand;
    
    idxXAxis = get(handles.pupDensXAxis, 'Value');
    idxYAxis = get(handles.pupDensYAxis, 'Value');
    
    set(handles.pupDensXAxis, 'Value', idxYAxis);
    set(handles.pupDensYAxis, 'Value', idxXAxis);

    pupDensAxis_Callback;
end

% --- Executes on selection change in lstIntGates.
function lstIntGates_Callback(~, ~, ~)
lstGates_Callback;
end

% --- Executes on selection change in pupKNN.
function pupKNN_Callback(hObject, eventdata, handles)
    plotScatter;
end

function txtRegexps_Callback(hObject, eventdata, handles)
    select_channels(get(hObject, 'String'));
end
% --- Executes on button press in btnAddRegexp.
function btnAddRegexp_Callback(hObject, eventdata, handles)
handles = gethand;
    regexps = retr('regexps');
    new_regexp = get(handles.txtRegexps, 'String');
    if isempty(new_regexp) || strcmp(new_regexp, '')
        % save the currently selected channels
        names = retr('channelNames');
        names = names';
        selected_channels = get(handles.lstChannels, 'Value');
        selected_names = names(selected_channels);
        selected_names = strrep(selected_names, '(', '\(');
        selected_names = strrep(selected_names, ')', '\)');
        selected_names = strcat(selected_names, '$');
        selected_names = strcat('^', selected_names);
        new_regexp = strjoin(selected_names, '|');
    end
    regexps{end+1, 1} = new_regexp;
    regexps{end, 2} = new_regexp;
	put('regexps', regexps);
    
    set(handles.lstRegexps, 'String', regexps(:, 1));
end

% --- Executes during object creation, after setting all properties.
function btnLoadGates_CreateFcn(hObject, ~, ~)
    placeIcon(hObject, 'add.gif');
end

% --- Executes during object creation, after setting all properties.
function btnSaveGate_CreateFcn(hObject, ~, ~)
    placeIcon(hObject, 'save.gif');
end

% --- Executes during object creation, after setting all properties.
function btnRemoveGate_CreateFcn(hObject, ~, ~)
    placeIcon(hObject, 'delete.gif');
end

% --- Executes during object creation, after setting all properties.
function btnReorderGateUp_CreateFcn(hObject, ~, ~)
    placeIcon(hObject, 'arrowup.gif');
end

% --- Executes during object creation, after setting all properties.
function btnHistIntLeft_CreateFcn(hObject, ~, ~)
    placeIcon(hObject, 'arrowleft_inac.png');
end

% --- Executes during object creation, after setting all properties.
function btnHistIntRight_CreateFcn(hObject, ~, ~)
    placeIcon(hObject, 'arrowright.png');
end

% --- Executes during object creation, after setting all properties.
function btnReorderGateDown_CreateFcn(hObject, ~, ~)
    placeIcon(hObject, 'arrowdown.gif');
end

% --- Executes during object creation, after setting all properties.
function btnDiff_CreateFcn(hObject, ~, ~)
    placeIcon(hObject, 'diff.gif');
end

% --- Executes during object creation, after setting all properties.
function btnSplitCellCycle_CreateFcn(hObject, ~, ~)
    placeIcon(hObject, 'cycle.png');
end

function btnIntersect_CreateFcn(hObject, ~, ~)
    placeIcon(hObject, 'intersect.png');
end

% --- Executes during object creation, after setting all properties.
function lstIntGates_CreateFcn(hObject, ~, ~)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end

% --- Executes during object creation, after setting all properties.
function pupKNN_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end

% --- Executes during object creation, after setting all properties.
function txtRegexps_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end

% --- Executes during object creation, after setting all properties.
function btnAddRegexp_CreateFcn(hObject, eventdata, handles)
    placeIcon(hObject, 'add.gif');
end

% --- Executes during object creation, after setting all properties.
function btnDensSwitchAxis_CreateFcn(hObject, eventdata, handles)
    placeIcon(hObject, 'switch1.gif');
end

function btnSwitchAxis_Callback(axis)
    handles = gethand;    
    
    idxYAxis = get(handles.pupYAxis, 'Value');

    if (axis==0) % == 'xy')
        idxOtherAxis = get(handles.pupXAxis, 'Value');
        set(handles.pupXAxis, 'Value', idxYAxis);
    elseif get(handles.pupZAxis, 'Value') > 1
        idxOtherAxis = get(handles.pupZAxis, 'Value')-1;
        set(handles.pupZAxis, 'Value', idxYAxis+1);
    end
    
    set(handles.pupYAxis, 'Value', idxOtherAxis);
    pupAxis_Callback;
end

function placeIcon(hObject, filename);
    try 
        [image_pic, map] = imread(filename);
        if ~isempty( map ) 
            image_pic = ind2rgb( image_pic, map ); 
        end 
        set(hObject,'cdata',image_pic); 
    catch e
        fprintf('warning: failed reading icon image ''%s''', filename);
    end

end

% --- Executes during object creation, after setting all properties.
function pupZAxis_CreateFcn(hObject, eventdata, handles)
    if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
        set(hObject,'BackgroundColor','white');
    end
end

function pupMarkerSize_CreateFcn(hObject, eventdata, handles)
    if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
        set(hObject,'BackgroundColor','white');
    end
end

% --- Executes on selection change in lstRegexps.
function lstRegexps_Callback(hObject, eventdata, handles)
handles = gethand;
sel_regexps = get(hObject, 'Value');
if (isempty(sel_regexps) || numel(sel_regexps) > 1)
    set(handles.lstChannels, 'Value', []);
else
    regexps = retr('regexps');
    new_reg = sprintf('(%s', regexps{sel_regexps(1), 2});
    for i=sel_regexps(2:end)
        new_reg = sprintf('%s|%s',new_reg,regexps{sel_regexps(i), 2});
    end
	new_reg = sprintf('%s)',new_reg);
    select_channels(new_reg);
    set(handles.txtRegexps, 'String', new_reg);
end
end

function select_channels(new_reg)
    handles = gethand;
    channel_names = retr('channelNames');
    if (isempty(channel_names)) return; end
        
    new_sel = find(cellfun(@(x)( ~isempty(x) ), regexp(channel_names, new_reg)));
    set(handles.lstChannels, 'Value', new_sel);
    lstChannels_Callback;
end

% --- Executes during object creation, after setting all properties.
function lstRegexps_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end

function cmiDeleteRegexp_callback
    handles = gethand;
    regexps = retr('regexps');

    regexps_names = get(handles.lstRegexps, 'String');
    sel_regexps = get(handles.lstRegexps, 'Value');
    
    regexps_names(sel_regexps) = [];
    regexps(sel_regexps, :) = [];
    
	set(handles.lstRegexps, 'Value', 1);
    set(handles.lstRegexps, 'String', regexps_names);
    put('regexps', regexps);
end

% --------------------------------------------------------------------
function cmnRegexps_Callback(hObject, eventdata, handles)
end

% --- Executes on button press in chkLegend.
function chkLegend_Callback(hObject, eventdata, handles)
    if ~get(hObject, 'Value')
        legend('off');
    else
        plotScatter;
    end
end

function sldrWanderlustWindowSize_Callback
    handles = gethand;
    ws = get(handles.sldrWanderlustWindowSize, 'Value');
    set(handles.txtWindowSize, 'String', sprintf('%1.2f', (ws-1)/100));
    plot_along_time;
end

% --------------------------------------------------------------------
function tSNE_each_gate_subsample(sample_size)
    handles = gethand;
    selected_gates = get(handles.lstGates, 'Value');

    for i=selected_gates
        set(handles.lstGates, 'Value', i);
        lstGates_Callback;
        runTSNE(2);
    end

end

function compare_tsne_subsample_runs

handles = gethand;
sessionData = retr('sessionData');
gateContext = retr('gateContext');
channel_names = retr('channelNames');
gates = retr('gates');
selected_gates = get(handles.lstGates, 'Value');
subsample_size = 6000;

% uniquely subsample 3 times the size of subsample_size/2
sample1 = randsample(gateContext, subsample_size/2);
sample2 = randsample(setdiff(gateContext, sample1), subsample_size/2);
sample_shared = randsample(setdiff(gateContext, union(sample1, sample2)), subsample_size/2);
createNewGate(sample1       , channel_names, {'sample1'});
createNewGate(sample2       , channel_names, {'sample2'});
createNewGate(sample_shared , channel_names, {'sample_shared'});

% update GUI with new gates 
gates = retr('gates');
% set(handles.lstGates, 'String', gates(:, 1));
% set(handles.lstIntGates, 'String', gates(:, 1));

% Add channel 'gate source'.
gate_source = zeros(size(sessionData,1), 1);

for j=selected_gates
    gate_source(gates{j, 2}) = j;
end
put('sessionData', sessionData);
addChannels({'gate_source'}, gate_source(:), 1:numel(gate_source), size(gates, 1)-2:size(gates, 1));

% run tsne twice on two pairs subsampled above.
nGates = size(gates, 1);
set(handles.lstGates, 'Value', [nGates-2 nGates]);
lstGates_Callback;
runTSNE(2);
set(handles.lstGates, 'Value', [nGates-1 nGates]);
lstGates_Callback;
runTSNE(2);

sessionData = retr('sessionData');
gates       = retr('gates');
ttsne1 = sessionData(:, end-3:end-2);
ttsne2 = sessionData(:, end-1:end);

sample_shared = gates{end, 2};

% compute pdist on the 2D tsne channels per source community
gate_sources = unique(sessionData(sample_shared, end-4))';
gate_source_corr = zeros(numel(gate_sources));
gate_source_count = zeros(numel(gate_sources));
for gate_source_idx = gate_sources
    gate_source_indices = find(sessionData(sample_shared, end-4) == gate_source_idx);
    gate_source_count( gate_source_idx ) = length( gate_source_indices );
    if( gate_source_count( gate_source_idx ) < 30 )
        continue;
    end
    
    gate_source_tsne1_pdist = pdist( ttsne1( sample_shared(gate_source_indices), : ) );
    gate_source_tsne2_pdist = pdist( ttsne2( sample_shared(gate_source_indices), : ) );
    
    % normalize the distances
    gate_source_tsne1_pdist = gate_source_tsne1_pdist./max(gate_source_tsne1_pdist);
    gate_source_tsne2_pdist = gate_source_tsne2_pdist./max(gate_source_tsne2_pdist);
    
    gate_source_corr( gate_source_idx ) = corr( gate_source_tsne1_pdist', gate_source_tsne2_pdist' );
end
gate_source_corr( gate_source_corr == 0 ) = NaN;

% plot the community size vs the pdist corrolation
current_figure = gcf;
h = figure;
hold on;
colors = distinguishable_colors(numel(gate_sources));
gate_names = remove_repeating_strings(gates(gate_sources, 1));
for gate_source_idx = gate_sources
    plot( gate_source_count(gate_source_idx), gate_source_corr(gate_source_idx), 'o', 'Color', colors(gate_source_idx,:) );
    text( gate_source_count(gate_source_idx), gate_source_corr(gate_source_idx), gate_names{gate_source_idx});
end
axis([0 1000 0 1]);
box on;
xlabel( 'size of manually gated cell type' );
ylabel( 'correlation (pearson)' );
% legend(gate_names);
n = now;
print(h, '-dpng', ['corr_plot_' date num2str(n) '.png']);
save(['corr_data_' date num2str(n)], 'sessionData', 'gates', 'gate_source_count', 'gate_source_corr');
hold off;

% Return to base figure
set(0,'CurrentFigure', current_figure);
end

function merge_gates
    handles = gethand;
    selected_gates = get(handles.lstGates, 'Value');
    sessionData = retr('sessionData');

    createNewGate(retr('gateContext'), retr('channelNames'));
    
    gates         = retr('gates');
        
	% if more than one gate was selected - add a channel 'gate source'
    if numel(selected_gates) > 1
        sessionData = retr('sessionData');
        v = zeros(size(sessionData,1), 1);

        for j=selected_gates
            v(gates{j, 2}) = j;
        end
        put('sessionData', sessionData);
        addChannels({'gate_source'}, v, 1:numel(v), size(gates, 1));
    end

end

function diff_gates
    handles = gethand;
    selected_gates = get(handles.lstGates, 'Value');
    sessionData = retr('sessionData');
    
    gates         = retr('gates');
    
    createNewGate(setdiff(gates{selected_gates(1), 2}, gates{selected_gates(2), 2}), retr('channelNames'));

%     set(handles.lstGates, 'String', gates(:, 1));
%     set(handles.lstIntGates, 'String', gates(:, 1));
    
end

function modularity_score=calculate_modularity
handles = gethand;
selected_gates = get(handles.lstGates, 'Value');
selected_channels = get(handles.lstChannels, 'Value');
sessionData = retr('sessionData');
gates = retr('gates');
gateContext	= retr('gateContext');

k_values = get(handles.pupKNN, 'String');
selection    = get(handles.pupKNN, 'Value');
k_neighbors = str2num(k_values{selection});
if (k_neighbors)
    str_k_neighbors = k_values{selection};
else
    str_k_neighbors = '5';
end

str_k_neighbors = inputdlg(...
            'enter K: ',...
            'Select K Neighbors', 1, {str_k_neighbors});
if (isempty(str_k_neighbors) || str2num(str_k_neighbors{1}) == 0)
    return;
end
k_neighbors = str2num(str_k_neighbors{1});
if numel(selected_gates) == 1
    [grouping_col,OK] = listdlg('PromptString','Select grouping channel: ',...
                                'Name', 'Grouping',...
                    'SelectionMode','single',...
                    'ListString',retr('channelNames'));

    if ~OK
        return;
    end
    grouping = sessionData(gateContext, grouping_col);
else
    grouping_col = 0;
    grouping = zeros(numel(gateContext));
    for i=selected_gates
        for j=1:numel(gateContext)
            if find(gates{i, 2} == gateContext(j))
                grouping(j) = i;
            end
        end
    end
end

disp 'Begin modularity computation ...';
mat = sessionData(gateContext, selected_channels);
disp 'Create KNN sparse matrix ...';
knn_sparse = create_k_sparse(mat, k_neighbors);

disp 'Calculate modularity ...';
modularity_score = spdists_modularity(knn_sparse, grouping)

disp 'Plotting ...';

set(handles.pupColorBy, 'Value', grouping_col+1);
set(handles.pupKNN, 'Value', k_neighbors+1);

current_figure = gcf;
h = figure;
plotScatter;
title (sprintf('Modularity: %g (KNN: %g)', modularity_score, k_neighbors));
print(h,'-dpng', ['corr_plot_' date num2str(now) '.png']);

% Return to base figure
set(0,'CurrentFigure', current_figure);

set(handles.pupKNN, 'Value', 1);

end

function txtWindowSize_Callback
    handles = gethand;
    window_size = get(handles.txtWindowSize, 'String');
    window_size = str2double(window_size);
    if ~isnan(window_size)
        window_size = (window_size*100)+1;
        sldrmin = get(handles.sldrWanderlustWindowSize, 'Min');
        sldrmax = get(handles.sldrWanderlustWindowSize, 'Max');
        window_size = in_range(sldrmin, sldrmax, window_size);
        set(handles.sldrWanderlustWindowSize, 'Value', window_size);
        plot_along_time;
    end
end

function num=in_range(min_r, max_r, num);
num = max(min_r, num);
num = min(max_r, num);
end
% --- Executes on slider movement.
function sldrBox_Callback(hObject, eventdata, handles)
    plotScatter;
end

% --- Executes during object creation, after setting all properties.
function sldrBox_CreateFcn(hObject, eventdata, handles)
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end
end

% --- Executes on slider movement.
function sldrThreshold_Callback(hObject, eventdata, handles)
    plotScatter
end

% --- Executes during object creation, after setting all properties.
function sldrThreshold_CreateFcn(hObject, eventdata, handles)
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end
end


% --- Executes on slider movement.
function sldrDiffKNN_Callback(hObject, eventdata, handles)
plotScatter;
end

% --- Executes during object creation, after setting all properties.
function sldrDiffKNN_CreateFcn(hObject, eventdata, handles)
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end
end

function txtDiffKNNThresh_Callback(hObject, eventdata, handles)
plotScatter;
end

% --- Executes during object creation, after setting all properties.
function txtDiffKNNThresh_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end

% --- Executes on slider movement.
function sldrDiffVoronoi_Callback(hObject, eventdata, handles)
plotScatter;
end

% --- Executes during object creation, after setting all properties.
function sldrDiffVoronoi_CreateFcn(hObject, eventdata, handles)
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end
end


% --- Executes during object creation, after setting all properties.
function rbgDiff_CreateFcn(hObject, eventdata, handles)
set(hObject, 'SelectionChangeFcn', @rbgDiff_Callback);
end

function rbgDiff_Callback(~, newSelection, oldSelection)
    handles = gethand;

    if (get(handles.rdbDiffBox, 'Value'))
        set(handles.pnlBoxDiffControls, 'Visible', 'on');
        set(handles.pnlKNNDiffControls, 'Visible', 'off');
        set(handles.pnlVoronoiDiffControls, 'Visible', 'off');
    elseif (get(handles.rdbDiffKNN, 'Value'))
        set(handles.pnlBoxDiffControls, 'Visible', 'off');
        set(handles.pnlKNNDiffControls, 'Visible', 'on');
        set(handles.pnlVoronoiDiffControls, 'Visible', 'off');
    else
        set(handles.pnlBoxDiffControls, 'Visible', 'off');
        set(handles.pnlKNNDiffControls, 'Visible', 'off');
        set(handles.pnlVoronoiDiffControls, 'Visible', 'on');
    end
    plotScatter;
end

function runPCA
    handles = gethand; % GUI handles to retrieve info from gui as below

    selected_gates    = get(handles.lstGates, 'Value'); % currently selected 
    selected_channels = get(handles.lstChannels, 'Value'); % ---------------

    session_data = retr('sessionData'); % all data
    gates       = retr('gates'); % all gates (names\indices) in cell array
    gate_context	= retr('gateContext'); % indices currently selected
    channel_names = retr('channelNames');
    
    coeff = pca(session_data(gate_context, selected_channels),2);
    addChannels({'PC1','PC2'}, coeff, gate_context); 
end

function add_gate(name, inds, channel_names)
gates = retr('gates');
gates{end+1, 1} = name;
gates{end, 2} = inds;
gates{end, 3} = channel_names;
put('gates', gates);
end

function flocked_data = flock(data)
    flocked_data = data;
    distance = 'mahalanobis';
    k = 15; 
	knn = parfor_spdists_knngraph( data, k, 'distance', distance, 'chunk_size', 1000, 'verbose', true );
	[li, lj, s] = find(knn);
    for i=1:size(data, 1)
        flocked_data(i, :) = median(data(li(lj==i), :));
    end

end
function create_cluster_means(opt_normalizemeans)
    handles = gethand; % GUI handles to retrieve info from gui as below

	selected_gates    = get(handles.lstGates, 'Value'); % currently selected
    selected_channel = get(handles.lstChannels, 'Value'); % ---------------

    session_data  = retr('sessionData'); % all data
    gates         = retr('gates'); 
    gate_context  = retr('gateContext'); % indices currently selected
    channel_names = retr('channelNames');

	data = session_data(gate_context, :);
    grouping = session_data(gate_context, selected_channel);
    
    if ~exist('opt_normalizemeans','var')
        selection = questdlg('Use normalized means?' ,...
                         'Normalize',...
                         'Yes','No','Yes');

        normalizemeans = strcmp(selection,'Yes');
    else
        normalizemeans = opt_normalizemeans;
    end
    
    if normalizemeans
        data = mynormalize(data, 99);
        data(:, selected_channel) = grouping;
    end
    
	cluster_ids = unique(grouping)';
    centers = zeros(length(cluster_ids), size(session_data, 2));
    i=0;
    for c=cluster_ids
        i = i+1;
        cluster_inds = find(grouping==c);
        centers(i, :) = mean(data(cluster_inds, :));
    end

    samplePercent = countmember(cluster_ids, grouping);
    samplePercent = samplePercent./max(samplePercent);
      
    session_data(end+1:end+size(centers,1), :) = centers;
    put('sessionData', session_data);
    
    
    cluster_gate_inds = size(session_data, 1)-(size(centers,1)-1):size(session_data, 1);
    
    if (numel(selected_gates) == 1)
        add_gate(sprintf('%s clusters - %s',...
            gates{selected_gates, 1},...
            channel_names{selected_channel}),...
            cluster_gate_inds,...
            channel_names);
    else
        add_gate(sprintf('clusters - %s',...
            channel_names{selected_channel}),...
            cluster_gate_inds,...
            channel_names);
    end
    
    addChannels({'sample percent', 'sample size'},...
        [countmember(cluster_ids, grouping)' samplePercent'],...
        cluster_gate_inds,...
        size(gates, 1)+1);

end

function create_indicator_channels_from_gates

    handles = gethand; % GUI handles to retrieve info from gui as below

    session_data = retr('sessionData'); % all data
    gates       = retr('gates'); % all gates (names\indices) in cell array

    selected_gates    = get(handles.lstGates, 'Value'); % currently selected 
    selected_gate_names = remove_repeating_strings(gates(selected_gates, 1));

    % generate descrete indicator channels for seelcted gates
    new_names = {};
    new_channels = zeros(size(session_data, 1), numel(selected_gates));
    for i=1:numel(selected_gates)
        new_names{end+1} = sprintf('%s - indicator', selected_gate_names{i});     
        new_channels(gates{selected_gates(i), 2}, i) =  1;
    end
    addChannels(new_names, new_channels, 1:size(session_data, 1));  
end

function louvainEach
    handles = gethand; % GUI handles to retrieve info from gui as below

    selected_gates    = get(handles.lstGates, 'Value'); % currently selected 
    selected_channels = get(handles.lstChannels, 'Value'); % ---------------

    session_data  = retr('sessionData'); % all data
    gates         = retr('gates'); % all gates (names\indices) in cell array
    gate_context  = retr('gateContext'); % indices currently selected
    channel_names = retr('channelNames');
    
    nKNeighbors = 15;
    
    for i=1:numel(selected_gates)
        % create a sparse KNN matrix
        sparse_adjacency_matrix = spdists_knngraph(session_data(gates{selected_gates(i), 2}, selected_channels), nKNeighbors, 'cosine', 5000);

        % Call louvain's matlab implementation
        [cmty mod] = spdists_louvain(sparse_adjacency_matrix);
        
        louvain_res = {cmty, mod};
        numel(mod)
        numel(cmty)
        
        % add new louvain channels
        mods = louvain_res{2};

        new_channel_names = cell(1, min (numel(mods), numel(louvain_res{1})));
        for ni=1:numel(new_channel_names)
            new_channel_names{ni} = sprintf('louvain_%g', mods(ni));
        end

        addChannels(new_channel_names, cell2mat(louvain_res{1}), gates{selected_gates(i), 2}, selected_gates(i));
    end
    
    return;

end

function phenoEach
    handles = gethand; % GUI handles to retrieve info from gui as below

    selected_gates    = get(handles.lstGates, 'Value'); % currently selected 
    selected_channels = get(handles.lstChannels, 'Value'); % ---------------

    session_data  = retr('sessionData'); % all data
    gates         = retr('gates'); % all gates (names\indices) in cell array
    gate_context  = retr('gateContext'); % indices currently selected
    channel_names = retr('channelNames');
 
    % ask user for cluster count to recognize
    k_neigh = inputdlg('Enter number of neigh: ',...
                          'cluster each', 1, {num2str(floor(sqrt(numel(gate_context)/2)))}); 
    
    if (isempty(k_neigh) || str2num(k_neigh{1}) == 0)
        return;
    end
    
    k_neigh = str2num(k_neigh{1});
    
    for i=1:numel(selected_gates)
        data = session_data(gates{selected_gates(i), 2}, selected_channels);
        
        [IDX, ~] = phenograph(data, k_neigh);
        addChannels({'pheno'}, IDX, gates{selected_gates(i), 2}, selected_gates(i));
    end
    
    return;  

end

function kmeansEach
    handles = gethand; % GUI handles to retrieve info from gui as below

    selected_gates    = get(handles.lstGates, 'Value'); % currently selected 
    selected_channels = get(handles.lstChannels, 'Value'); % ---------------

    session_data  = retr('sessionData'); % all data
    gates         = retr('gates'); % all gates (names\indices) in cell array
    gate_context  = retr('gateContext'); % indices currently selected
    channel_names = retr('channelNames');
 
    for i=1:numel(selected_gates)
        k_clusters = 4;
        addChannels({'kmeans4'}, kmeans(session_data(gates{selected_gates(i), 2}, selected_channels), k_clusters,'Display', 'iter', 'EmptyAction', 'singleton'), gates{selected_gates(i), 2}, selected_gates(i));
    end
    
    return;  

end

% ------------
% compare between two maps
% ------------
function compareMaps

    handles = gethand; % GUI handles to retrieve info from gui as below
    session_data  = retr('sessionData'); % all data
    gate_context  = retr('gateContext'); % indices currently selected
    channel_names = retr('channelNames');

    map_compareGUI('session_data',session_data(gate_context,:),'channelNames',channel_names)
end

% ------------
% add your own code here. please refer to examples below.
% ------------
function openEndedAction

    handles = gethand; % GUI handles to retrieve info from gui as below

    selected_gates    = get(handles.lstGates, 'Value'); % currently selected 
    selected_channels = get(handles.lstChannels, 'Value'); % ---------------

    session_data  = retr('sessionData'); % all data
    gates         = retr('gates'); % all gates (names\indices) in cell array
    gate_context  = retr('gateContext'); % indices currently selected
    channel_names = retr('channelNames');
        
    
     return;
    
   %% Reorder thymus panel channel names
    chDNA = cellstrfnd(channel_names, 'dna');
    channel_names(chDNA) = strcat('..', channel_names(chDNA))

    chBEAD = cellstrfnd(channel_names, 'bead');
    channel_names(chBEAD) = strcat('.', channel_names(chBEAD))

    chBEAD = cellstrfnd(channel_names, 'live');
    channel_names(chBEAD) = strcat('..', channel_names(chBEAD))

    channel_names{1} =  '...time';
    channel_names{2} =  '...viability';
    
    [idx, sorted] = sort(channel_names);
    session_data = session_data(:, sorted);
    for g=1:size(gates, 1)
        gates{g, 3} =   idx;
    end
    
    put('gates', gates);
    put('sessionData', session_data);

    data = session_data(gate_context, selected_channels);
    ell = GenerateLoop(data(:, 1)', data(:,2)');
    save 'ell3.mat' ell;
    IDX = knnsearch(ell, data);
    
    addChannels({'ERA'}, IDX);

%     return;
    
    %% gate on DNA vs DNA channels
    cofactor = 5;
    gate_ind = selected_gates(1);
    debris_threshold = 0.9;

    chDNA = cellstrfnd(gates{gate_ind, 3}, 'DNA');
    DNA = session_data(gates{gate_ind, 2}, chDNA);
    DNA = asinh( DNA ./ cofactor ); % transform data
    dna_gm = gmdistribution.fit(DNA, 2);
    [~,j] = max(sum(dna_gm.mu,2));
    P = dna_gm.posterior(DNA);
    debris = gates{gate_ind, 2}(P(:,j) > debris_threshold);
    addGate('Thymus1 removed debris by dna', debris);
%     return; 
    
    %% create a hsin transformed gate of the data
    gates = retr('gates');
	selected_gates    = get(handles.lstGates, 'Value'); % currently selected 
    cofactor = 5;
    insert_loc = size(session_data, 1)+1;
    session_data(end+1:end+length(gate_context), :) =  asinh( session_data(gate_context, :) ./ cofactor );
    put('sessionData', session_data);
    addGate([gates{selected_gates, 1} ' transformed'], insert_loc:(insert_loc+length(gate_context)-1));
%     return;

    %% bubble check
    binsize = 20;
    current_figure = gcf;
    h = figure;
    try
        for i=1:numel(selected_gates)
            data = session_data(gates{selected_gates(i), 2}, :);
            nbins = floor(size(data,1)/binsize);
            data = data(1:nbins*binsize, :);
            [~, idx] = sort(data(:, 1));
            for ch=2:length(channel_names)
                valuesperbin = reshape(data(idx, ch), binsize, nbins);
                means = mean(valuesperbin, 1);
                plot(1:nbins, means , '-');
                title(channel_names{ch});
                print(h,'-dpng', validfilename(sprintf('%s %s vs time', gates{i, 1}, channel_names{ch})));
            end
        end
    catch e

        % Return to base figure
        set(0,'CurrentFigure', current_figure);
        close(h);

        disp(getReport(e,'extended'));
    end
    set(0,'CurrentFigure', current_figure);
    close(h);

%     return;

   %% grab an untrasformed gate from transformed gate
    transformed_cd8_gate_indices = gates{selected_gates(2), 2};
    transformed_all_points = gates{selected_gates(1), 2};
    gate_inds = find(ismember(transformed_all_points,transformed_cd8_gate_indices));
    CD8high_untransformed = gates{selected_gates(1) - 1, 2}(gate_inds);
    addGate('cd8high untransformed', CD8high_untransformed);
    return;
 
    %%
%     data = session_data(gate_context, selected_channels);
%     IDX = knnsearch(ell, data);
%     addChannels({'ERA_all'}, IDX); 
%     IDX = knnsearch(ell(1:65, :), data);
%     addChannels({'ERA_all'}, IDX); 
%     return;
    data = session_data(gate_context, selected_channels);
    ell = GenerateLoop(data(:, 1)', data(:,2)');
    save 'ell.mat' ell;
    IDX = knnsearch(ell, data);
    
    addChannels({'ERA'}, IDX);

    return;
    
    axis tight;
    OptionZ.FrameRate=30;OptionZ.Duration=8;OptionZ.Periodic=true;
    CaptureFigVid([-130,15;-165,10;-200,80;-235,15;-280,10],'3DAML',OptionZ)

    
    highprc = prctile(session_data(gate_context, selected_channels), 90);
    channel_names{selected_channels(highprc > 1.5)}
    return;
    phenoEach;
    return;
    
    for gi=1:numel(selected_gates)
        set(handles.lstGates, 'Value', selected_gates(gi)); % currently selected 
        lstGates_Callback;
        create_cluster_means(false);
    end
    return;

	create_indicator_channels_from_gates
    return;

    % create a knn graph and export cluster means to graphml
    [s,v] = listdlg('PromptString','Select a cluster channel:',...
            'SelectionMode','single',...
            'InitialValue', length(channel_names),...
            'ListString',channel_names);
    if ~v
        return;
    end
    cluster_channel = s;
    
    data = session_data(gate_context, selected_channels);
    grouping = session_data(gate_context, cluster_channel);
    
    selection = questdlg('Use normalized means?' ,...
                     'Normalize',...
                     'Yes','No','Yes');
    if strcmp(selection,'Yes')
        data = mynormalize(data, 99);
    end

    cluster_ids = unique(grouping)';
    centers = zeros(length(cluster_ids), length(selected_channels));
    i=0;
    for c=cluster_ids
        i = i+1;
        cluster_inds = find(grouping==c);
        centers(i, :) = mean(data(cluster_inds, :));
    end
    
    dists         = pdist2(centers, centers, 'euclidean','Smallest',3);
    dists         = 1 - dists./max(max(dists));
    
    sampleID      = cluster_ids;
    samplePercent = countmember(cluster_ids, grouping);
    samplePercent = samplePercent./max(samplePercent);
    means         = centers;
    names         = channel_names(selected_channels);
    
    node_names  = cell(1, length(cluster_ids));
    for i = 1:length(cluster_ids)
        for ch=1:length(selected_channels)
            node_names{i} = sprintf('%s %s=%g ', node_names{i}, channel_names{selected_channels(ch)}, centers(i, ch));
        end
    end
    
	coeff = pca(means,2);
	scatter_by_point(coeff(:,1), coeff(:,2), cluster_ids', floor(samplePercent*2000)); %plotting
    axis([-1 1 -1 1])

%   write_graphml(dists, sampleID, samplePercent, means, names, node_names);
    
    return;
    
    % create a knn graph and export cluster means to graphml
    [s,v] = listdlg('PromptString','Select a cluster channel:',...
            'SelectionMode','single',...
            'InitialValue', length(channel_names),...
            'ListString',channel_names);
    if ~v
        return;
    end
    cluster_channel = s;
    
    data = session_data(gate_context, selected_channels);
    grouping = session_data(gate_context, cluster_channel);
    
    selection = questdlg('Use normalized means?' ,...
                     'Normalize',...
                     'Yes','No','Yes');
    if strcmp(selection,'Yes')
        data = mynormalize(data, 99);
    end

    cluster_ids = unique(grouping)';
    centers = zeros(length(cluster_ids), length(selected_channels));
    i=0;
    for c=cluster_ids
        i = i+1;
        cluster_inds = find(grouping==c);
        centers(i, :) = mean(data(cluster_inds, :));
    end
    
    dists         = pdist2(centers, centers, 'euclidean','Smallest',3);
    dists         = 1 - dists./max(max(dists));
    
    sampleID      = cluster_ids;
    samplePercent = countmember(cluster_ids, grouping);
    samplePercent = samplePercent./max(samplePercent);
    means         = centers;
    names         = channel_names(selected_channels);
    
    node_names  = cell(1, length(cluster_ids));
    for i = 1:length(cluster_ids)
        for ch=1:length(selected_channels)
            node_names{i} = sprintf('%s %s=%g ', node_names{i}, channel_names{selected_channels(ch)}, centers(i, ch));
        end
    end
    
	coeff = pca(means,2);
	scatter_by_point(coeff(:,1), coeff(:,2), cluster_ids, floor(samplePercent*32)); %plotting

%     write_graphml(dists, sampleID, samplePercent, means, names, node_names);
    
    return;
    
    % ask user for cluster count to recognize
    n_clusters = inputdlg('Enter number of clusters: ',...
                          'EMGM', 1, {num2str(floor(sqrt(numel(gate_context)/2)))}); 
    
    if (isempty(n_clusters) || str2num(n_clusters{1}) == 0)
        return;
    end

    k_clusters = str2num(n_clusters{1});
    
    for gatei = selected_gates
        % run EMGM
        clust_alg = 'EMGM';
        [IDX, model, llh] = emgm(session_data(gates{gatei, 2}, selected_channels)', k_clusters);
        IDX = IDX';

        % add results to GUI
        addChannels({sprintf('%s%g',clust_alg, k_clusters)}, IDX, gates{gatei, 2}, gatei);
    end
    
    return;
    
    % generate descrete indicator channels from a grouping\cluster channel
    new_names = {};
    new_channels = zeros(size(session_data, 1), numel(unique(session_data(:, selected_channels))));
    for i=unique(session_data(:, selected_channels))'
        new_names{i+1} = sprintf('%g - phase', i);     
        new_channels(:, i+1) =  (session_data(:, selected_channels) == i);
    end
    addChannels(new_names, new_channels);
    return;
%     data = session_data(gate_context, selected_channels);
%     flocked_data = flock(data);
% 
%     new_channel_names = strcat(channel_names(selected_channels), {' flocked'});
%     addChannels(new_channel_names, flocked_data, gate_context);
%     return;
    
%     BW1 = edge(data,'sobel');
%     imshow(BW1);
%     return;

% 	[~, channel_names] = xlsread('/Users/mtadmor/Google Drive/Mitosis/Gabri/16-Feb-2014_CellCCycleTimecourse_Crowding_featurelist.xls');
%     channel_names = channel_names'; 
%     gates{1, 3} = channel_names;
%     put('gates', gates);
%     
%    % plotEnrichment;
%      return;
 
        % === run wonderlust
        
        normalize = true;
        data = session_data(gate_context, selected_channels);

        if (normalize)
            disp('Normalizing according to the 99th percentile...');
            data = data-repmat(prctile(data, 1, 1), size(data,1),1);
            data = data./repmat(prctile((data), 99, 1),size(data,1),1);

            data(data > 1) = 1;
            data(data < 0) = 0;
            data(isinf(data)) = 0;

        end
        distance = 'mahalanobis';
        k = 10; 
        l = 40;
        num_graphs = 1; 

    %     [x,y]=ginput
    %     some_channel = 59;
    %     [~, s] = max(session_data(gate_context, some_channel));
%         s = knnsearch(session_data(gates{end, 2}(1), selected_channels), session_data(gates{selected_gates(1),2}, selected_channels));
%         [~, s] = ismember(data, session_data(gates{end, 2}(1), selected_channels), 'rows');
        s = find(gates{end, 2}(1)==gate_context);

        num_landmarks = 30; % numel(gate_context)/(k^3);

        predefined_l = [];
%         timepoint_marker = 115;
%         for t=unique(session_data(gate_context, timepoint_marker))'
%             ts = find(session_data(gate_context, timepoint_marker) == t);
%             predefined_l(end+1) = randsample(ts, 1);
%         end
        [traj, lnn] = wanderlust(data, distance, k, l, num_graphs, s, num_landmarks, true, predefined_l);

        save('lnn', 'lnn');
        traj = traj';
 
        new_channel_names = strcat({'traj '},int2str((1:num_graphs).'));

        addChannels(new_channel_names, traj, gate_context);
   
    return;
     
%     % change gate names
%     new_names = remove_repeating_strings(gates(selected_gates, 1));
%     gates(selected_gates, 1) = cellfun(@(name)strrep(name, '-channel 120-', '-'), new_names, 'UniformOutput', false)
%     put('gates', gates);
%     
%     return;

mi = zeros(1, 111);
    for i=1:111
        try
    mi(i) = compute_dremi([T phi(:, i)], .9);
        catch 
            continue;
        end
    end
    
    N = length(gate_context);
    learn_inds = randsample(1:N,floor(.2*N));

    phi = session_data(gate_context, selected_channels);
    T = session_data(gate_context, 4); %+ noise(learn_inds,ind, 1);

    lambda = 5;
    iters = 100;

    flag1 = 1; % type
    flag2 = 0; % default gamma, init to find min. solution
    flag3 = 1; % verbose


    [mu,dmu,k,gamma] = sparse_learning(phi, T, lambda, iters, flag1, flag2, flag3);
    
    return;

    % --------- examples and explanation ----------------------------------
    % ---- for example, to plot the curently selected gates\channels, write 
    % along these lines:
    % plot(session_data(gateContext, selected_channels));
    %
    %
    % ---- to plot in a new window and save figure the expression of a list of
    % channels ([1 2 7 12 13]) you can write along these lines:
    %
    %         current_figure = gcf;
    %         h = figure;
    % 
    %         surface_markers = [1 2 7 12 13 15];
    %         channel_names = gates(selected_gates(1), 3);
    %         try
    %             for j=1:numel(surface_markers)
    %                 subplot(2, 3, j);
    %                 data_space = session_data(gate_context, selected_channels);
    %                 vX = data_space(:,1);
    %                 vY = data_space(:,2);
    %                 vColor = sessionData(gate_context, surface_markers(j));
    %                 scatter(vX, vY, 6, vColor, '.');
    % 
    %                 title(channel_names{surface_markers(j)}, 'FontSize', 14);
    %                 box on;
    %                 set(gca,'xtick',[],'ytick',[]);
    %                 axis([min(vX) max(vX) min(vY) max(vY)]);
    %             end
    %             set(gcf, 'PaperPositionMode', 'auto');
    %             print(h,'-dpng', gates{selected_gates(1), 1});
    % 
    %             % Return to base figure
    %             set(0,'CurrentFigure', current_figure);
    %             close(h);        
    %         catch err
    %             display(err);
    %             
    %             % Return to base figure
    %             set(0,'CurrentFigure', current_figure);
    %             close(h);
    %         end
    %
    %  --------- end of examples --------------------------------------------



    %     thetwelve = [5 7 9 10 11 12 14 17 18 21 28 30];
    %     printColorExpressionInSelectedGatesOverSelectedChannels(thetwelve);

    %      for selected_channel=selected_channels
    %         current_figure = gcf;
    %         h = figure;
    %         
    %         try
    %              for ind=1:numel(selected_gates)
    %                 subplot(2, 3, ind);
    %                 data = session_data(gates{selected_gates(ind),2}, selected_channel);
    %                 dplot(data);
    %                 box on;
    %                 set(gca,'ytick',[]);
    %                 axis tight;
    %                 healthies = [2 3 5 6 7];
    %                 title(sprintf('H%g', healthies(ind)));
    %              end
    %              print(h,'-dpng', channel_names{selected_channel});
    %         catch err
    %            
    %             % Return to base figure
    %             set(0,'CurrentFigure', current_figure);
    %             close(h);
    %             
    %             rethrow(err);
    %         end
    %         set(0,'CurrentFigure', current_figure);
    %         close(h);
    %      end
end

