function varargout = louvain(varargin)
% LOUVAIN MATLAB code for louvain.fig
%      LOUVAIN by itself, creates a new LOUVAIN or raises the
%      existing singleton*.
%
%      H = LOUVAIN returns the handle to a new LOUVAIN or the handle to
%      the existing singleton*.
%
%      LOUVAIN('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in LOUVAIN.M with the given input arguments.
%
%      LOUVAIN('Property','Value',...) creates a new LOUVAIN or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before louvain_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to louvain_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help louvain

% Last Modified by GUIDE v2.5 03-May-2012 15:59:46

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @louvain_OpeningFcn, ...
                   'gui_OutputFcn',  @louvain_OutputFcn, ...
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
% --- Executes just before louvain is made visible.
function louvain_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to louvain (see VARARGIN)

% Choose default command line output for louvain
handles.output = {};

% Update handles structure
guidata(hObject, handles);
setappdata(0,'hlouv',gcf);

% Insert custom Title and Text if specified by the user
% Hint: when choosing keywords, be sure they are not easily confused 
% with existing figure properties.  See the output of set(figure) for
% a list of figure properties.

hgui=getappdata(0,'hlouv');
setappdata(hgui, 'data', varargin{1});

if(nargin > 5)
    for index = 3:2:(nargin-3),
        if nargin-3==index, break, end
        switch lower(varargin{index})
         case 'title'
          set(hObject, 'Name', varargin{index+1});
        case 'channelnames'
          set(handles.lstChannels, 'String', varargin{index+1});
          set(handles.lstChannels, 'Value', varargin{2});
        end
    end
else
    handles.output = runLouvain;

    % Update handles structure
    guidata(hObject, handles);

    % Use UIRESUME instead of delete because the OutputFcn needs
    % to get the updated handles structure.
    uiresume(handles.figure1);
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

% UIWAIT makes louvain wait for user response (see UIRESUME)
uiwait(handles.figure1);
end
% --- Outputs from this function are returned to the command line.
function varargout = louvain_OutputFcn(hObject, eventdata, handles)
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
        handles.output = {1, [], 0,0};
        % The GUI is no longer waiting, just close it
        delete(hObject);
    end
end

% --- Executes on key press over figure1 with no controls selected.
function figure1_KeyPressFcn(hObject, eventdata, handles)

    % Check for "enter" or "escape"
    if isequal(get(hObject,'CurrentKey'),'escape')
        % User said no by hitting escape
        handles.output = {1, [], 0,0};

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

    [c, m] = runLouvain;
    handles.output = {c, m};

    % Update handles structure
    guidata(hObject, handles);

    % Use UIRESUME instead of delete because the OutputFcn needs
    % to get the updated handles structure.
    uiresume(handles.figure1);
end

function [cmty mod]=runLouvain
    hgui=getappdata(0,'hlouv');
    handles=guihandles(hgui);

    selectedChannels = get(handles.lstChannels, 'Value');
    nKNeighbors = str2num(get(handles.txtKNeighbors, 'String'));
    nMaxDistance = str2num(get(handles.txtMaxDistance, 'String'));
    hgui=getappdata(0,'hlouv');  
    data=getappdata(hgui, 'data');
    
    % create a sparse KNN matrix
    sparse_adjacency_matrix = spdists_knngraph(data(:, selectedChannels), nKNeighbors, 'euclidean', 5000);

    % Call louvain's matlab implementation
	[cmty mod] = spdists_louvain(sparse_adjacency_matrix);
            numel(mod)
        numel(cmty)

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