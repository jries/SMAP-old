function varargout = create_metaGUI(varargin)
% MATLAB code for create_meta.fig
%      
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Last Modified by GUIDE v2.5 03-May-2012 15:59:46

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @create_metaGUI_OpeningFcn, ...
                   'gui_OutputFcn',  @create_metaGUI_OutputFcn, ...
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
% --- Executes just before interface is made visible.
function create_metaGUI_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to LOUVAIN_meta_clusters (see VARARGIN)

% Choose default command line output for finding meta clusters
handles.output = {};

% Update handles structure
guidata(hObject, handles);
setappdata(0,'hwand',gcf);

% Insert custom Title and Text if specified by the user
% Hint: when choosing keywords, be sure they are not easily confused 
% with existing figure properties.  See the output of set(figure) for
% a list of figure properties.

hgui=getappdata(0,'hwand');

% Finding function inputs
for i=1:length(varargin)-1        
    
    if(strcmp(varargin{i}, 'title'))
        set(hObject, 'Name', varargin{i+1});
    elseif(strcmp(varargin{i}, 'params'))
        params = varargin{i+1};
        set(handles.txtNeighbors, 'String', num2str(params.neighbors));
        metrics = get(handles.pupMetric, 'String');
        metric = find(strcmpi(metrics, params.metric));
        set(handles.pupMetric, 'Value', metric);
        set(handles.lstGates, 'String', params.selected_gates);
        %set(handles.lstClusterChannel, 'Value', params.all_channels);
        set(handles.lstClusterChannel, 'String', params.all_channels);

%     elseif(strcmp(varargin{i}, 'gates'))
%         set(handles.lstGates, 'String', varargin{i+1});
    end
end

% Determine the position of the dialog - centered on the callback figure
% if available, else, centered on the screen
FigPos=get(0,'DefaultFigurePosition');
OldUnits = get(hObject, 'Units');
set(hObject, 'Units', 'pixels');
OldPos = get(hObject,'Position');
FigWidth = OldPos(3);
FigHeight = OldPos(4);
if isempty(gcbf)
    ScreenUnits=get(0,'Units');
    set(0,'Units','pixels');
    ScreenSize=get(0,'ScreenSize');
    set(0,'Units',ScreenUnits);

    FigPos(1)=1/2*(ScreenSize(3)-FigWidth);
    FigPos(2)=2/3*(ScreenSize(4)-FigHeight);
else
    GCBFOldUnits = get(gcbf,'Units');
    set(gcbf,'Units','pixels');
    GCBFPos = get(gcbf,'Position');
    set(gcbf,'Units',GCBFOldUnits);
    FigPos(1:2) = [(GCBFPos(1) + GCBFPos(3) / 2) - FigWidth / 2, ...
                   (GCBFPos(2) + GCBFPos(4) / 2) - FigHeight / 2];
end
FigPos(3:4)=[FigWidth FigHeight];
set(hObject, 'Position', FigPos);
set(hObject, 'Units', OldUnits);

% Show a question icon from dialogicons.mat - variables questIconData
% and questIconMap
load dialogicons.mat

IconData=questIconData;
questIconMap(256,:) = get(handles.figure1, 'Color');
IconCMap=questIconMap;

set(handles.figure1, 'Colormap', IconCMap);

% Make the GUI modal
set(handles.figure1,'WindowStyle','modal')

% UIWAIT makes meta cluster finding wait for user response (see UIRESUME)
uiwait(handles.figure1);
end
% --- Outputs from this function are returned to the command line.
function varargout = create_metaGUI_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;

% The figure can be deleted now
delete(handles.figure1);
end

% --- Executes when user attempts to close figure1.
function figure1_CloseRequestFcn(hObject, eventdata, handles)

    if isequal(get(hObject, 'waitstatus'), 'waiting')
        % The GUI is still in UIWAIT, us UIRESUME
        uiresume(hObject);
    else
        handles.output = [];
        % The GUI is no longer waiting, just close it
        delete(hObject);
    end
end

% --- Executes on key press over figure1 with no controls selected.
function figure1_KeyPressFcn(hObject, eventdata, handles)

    % Check for "enter" or "escape"
    if isequal(get(hObject,'CurrentKey'),'escape')
        % User said no by hitting escape
        handles.output = [];

        % Update handles structure
        guidata(hObject, handles);

        uiresume(handles.figure1);
    end    

    if isequal(get(hObject,'CurrentKey'),'return')
        uiresume(handles.figure1);
    end   
end

% --- Executes on button press in btnRun.
function btnRun_Callback(hObject, eventdata, handles)

    W = findMetaClusters;
        
    handles.output = W;

    % Update handles structure
    guidata(hObject, handles);

    % Use UIRESUME instead of delete because the OutputFcn needs
    % to get the updated handles structure.
    uiresume(handles.figure1);
end

function W=findMetaClusters
    hgui=getappdata(0,'hwand');
    handles=guihandles(hgui);

    W.selected_gates = get(handles.lstGates, 'Value');
    W.cluster_channel = get(handles.lstClusterChannel, 'Value');
    W.neighbors = str2num(get(handles.txtNeighbors, 'String'));

    nMetric = get(handles.pupMetric, 'Value');
    strMetrics = get(handles.pupMetric, 'String');
    W.metric = strMetrics{nMetric};
    W.all_channels = get(handles.lstClusterChannel, 'String');
            
    
end

% --- Executes on button press in btnCancel.
function btnCancel_Callback(hObject, eventdata, handles)
% hObject    handle to btnCancel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles.output = {};

% Update handles structure
guidata(hObject, handles);

% Use UIRESUME instead of delete because the OutputFcn needs
% to get the updated handles structure.
uiresume(handles.figure1);
end

function txtMaxDistance_Callback(hObject, eventdata, handles)
% hObject    handle to txtMaxDistance (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of txtMaxDistance as text
%        str2double(get(hObject,'String')) returns contents of txtMaxDistance as a double
end

% --- Executes during object creation, after setting all properties.
function txtMaxDistance_CreateFcn(hObject, eventdata, handles)
% hObject    handle to txtMaxDistance (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end


function txtKNeighbors_Callback(hObject, eventdata, handles)
% hObject    handle to txtKNeighbors (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of txtKNeighbors as text
%        str2double(get(hObject,'String')) returns contents of txtKNeighbors as a double
end

% --- Executes during object creation, after setting all properties.
function txtKNeighbors_CreateFcn(hObject, eventdata, handles)
% hObject    handle to txtKNeighbors (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end

% --- Executes on selection change in lstChannels.
function lstChannels_Callback(hObject, eventdata, handles)
% hObject    handle to lstChannels (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns lstChannels contents as cell array
%        contents{get(hObject,'Value')} returns selected item from lstChannels
end

% --- Executes during object creation, after setting all properties.
function lstChannels_CreateFcn(hObject, eventdata, handles)
% hObject    handle to lstChannels (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end