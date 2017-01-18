function varargout = addROIGUI(varargin)
% addROIGUI M-file for addROIGUI.fig
%      addROIGUI, by itself, creates a new addROIGUI or raises the existing
%      singleton*.
%
%      H = addROIGUI returns the handle to a new addROIGUI or the handle to
%      the existing singleton*.
%
%      addROIGUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in addROIGUI.M with the given input arguments.
%
%      addROIGUI('Property','Value',...) creates a new addROIGUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before addROIGUI_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to addROIGUI_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES
%
% Copyright (C) 2016, Danuser Lab - UTSouthwestern 
%
% This file is part of u-track.
% 
% u-track is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
% 
% u-track is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
% 
% You should have received a copy of the GNU General Public License
% along with u-track.  If not, see <http://www.gnu.org/licenses/>.
% 
% 

% Edit the above text to modify the response to help addROIGUI

% Last Modified by GUIDE v2.5 29-Apr-2015 14:56:21

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @addROIGUI_OpeningFcn, ...
                   'gui_OutputFcn',  @addROIGUI_OutputFcn, ...
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


% --- Executes just before addROIGUI is made visible.
function addROIGUI_OpeningFcn(hObject,eventdata,handles,varargin)

% Check input
% The mainFig and procID should always be present
% procCOnstr and procName should only be present if the concrete process
% initation is delegated from an abstract class. Else the constructor will
% be directly read from the package constructor list.
ip = inputParser;
ip.addRequired('hObject',@ishandle);
ip.addRequired('eventdata',@(x) isstruct(x) || isempty(x));
ip.addRequired('handles',@isstruct);
ip.addOptional('MD',[],@(x)isa(x,'MovieData'));
ip.addParamValue('mainFig',-1,@ishandle);
ip.parse(hObject,eventdata,handles,varargin{:});

userData.MD =ip.Results.MD;
userData.mainFig =ip.Results.mainFig;

% Set up copyright statement
set(handles.text_copyright, 'String', getLCCBCopyright());

% Set up available input channels
set(handles.listbox_selectedChannels,'String',userData.MD.getChannelPaths(), ...
    'UserData',1:numel(userData.MD.channels_));

% Save the image directories and names (for cropping preview)
userData.nFrames = userData.MD.nFrames_;
userData.imPolyHandle.isvalid=0;
userData.ROI = [];
userData.previewFig=-1;
userData.helpFig=-1;

% Read the first image and update the sliders max value and steps
userData.chanIndex = 1;
set(handles.edit_frameNumber,'String',1);
set(handles.slider_frameNumber,'Min',1,'Value',1,'Max',userData.nFrames,...
    'SliderStep',[1/max(1,double(userData.nFrames-1))  10/max(1,double(userData.nFrames-1))]);
if userData.nFrames == 1
    set(handles.edit_frameNumber, 'Enable', 'off');
    set(handles.slider_frameNumber, 'Enable', 'off');
end
userData.imIndx=1;
userData.imData=mat2gray(userData.MD.channels_(userData.chanIndex).loadImage(userData.imIndx));
    
set(handles.listbox_selectedChannels,'Callback',@(h,event) update_data(h,event,guidata(h)));

userData_main = get(ip.Results.mainFig, 'UserData');
set(handles.figure1,'CurrentAxes',handles.axes_help);
Img = image(userData_main.questIconData);
set(gca, 'XLim',get(Img,'XData'),'YLim',get(Img,'YData'),...
    'visible','off','YDir','reverse');
set(Img,'ButtonDownFcn',@icon_ButtonDownFcn,...
    'UserData', struct('class', mfilename))

% Choose default command line output for addROIGUI
handles.output = hObject;

% Update user data and GUI data
set(hObject, 'UserData', userData);
guidata(hObject, handles);

% --- Outputs from this function are returned to the command line.
function varargout = addROIGUI_OutputFcn(~, ~, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;

% --- Executes on button press in pushbutton_cancel.
function pushbutton_cancel_Callback(~, ~, handles)
% Delete figure
delete(handles.figure1);

% --- Executes during object deletion, before destroying properties.
function figure1_DeleteFcn(hObject, ~, handles)
% Notify the package GUI that the setting panel is closed
userData = get(handles.figure1, 'UserData');

if ~isempty(userData)
    
    if ishandle(userData.helpFig), delete(userData.helpFig); end
    if ishandle(userData.previewFig), delete(userData.previewFig); end
    
    set(handles.figure1, 'UserData', userData);
    guidata(hObject,handles);
end

% --- Executes on key press with focus on pushbutton_save and none of its controls.
function pushbutton_save_KeyPressFcn(~, eventdata, handles)

if strcmp(eventdata.Key, 'return')
    pushbutton_done_Callback(handles.pushbutton_save, [], handles);
end

% --- Executes on key press with focus on figure1 and none of its controls.
function figure1_KeyPressFcn(~, eventdata, handles)

if strcmp(eventdata.Key, 'return')
    pushbutton_done_Callback(handles.pushbutton_save, [], handles);
end

% --- Executes on button press in checkbox_preview.
function update_data(hObject, eventdata, handles)
userData = get(handles.figure1, 'UserData');

% Retrieve the channel index
props=get(handles.listbox_selectedChannels,{'UserData','Value'});
chanIndex = props{1}(props{2});
imIndx = get(handles.slider_frameNumber,'Value');

% Load a new image if either the image number or the channel has been changed
if (chanIndex~=userData.chanIndex) ||  (imIndx~=userData.imIndx)
    % Update image flag and dat
    userData.imData=mat2gray(userData.MD.channels_(chanIndex).loadImage(imIndx));
    userData.updateImage=1;
    userData.chanIndex=chanIndex;
    userData.imIndx=imIndx;
    
    % Update roi
    if userData.imPolyHandle.isvalid
        userData.ROI=getPosition(userData.imPolyHandle);
    end
else
    userData.updateImage=0;
end


% Create figure if non-existing or closed
if ~isfield(userData, 'previewFig') || ~ishandle(userData.previewFig)
    userData.previewFig = figure('NumberTitle','off','Name',...
        'Select the region of interest','DeleteFcn',@close_previewFig,...
        'UserData',handles.figure1);
    axes('Position',[.05 .05 .9 .9]);
    userData.newFigure = 1;
else
    figure(userData.previewFig);
    userData.newFigure = 0;
end

% Retrieve the image object handle
imHandle =findobj(userData.previewFig,'Type','image');
if userData.newFigure || userData.updateImage
    if isempty(imHandle)
        imHandle=imshow(userData.imData);
        axis off;
    else
        set(imHandle,'CData',userData.imData);
    end
end

% set the userdata before interactive impoly in case the figure changes
% first
set(handles.figure1, 'UserData', userData);
guidata(hObject,handles);

if userData.imPolyHandle.isvalid
    % Update the imPoly position
    setPosition(userData.imPolyHandle,userData.ROI)
else
    % Do not create a new impoly if one is in progress
    if(isempty(findobj(gca,'Tag','impoly')))
    
        % Create a new imPoly object and store the handle
        % since impoly blocks, create the constraint function first in case
        % something happens later
        fcn = makeConstrainToRectFcn('impoly',get(imHandle,'XData'),get(imHandle,'YData'));
        if ~isempty(userData.ROI)
            userData.imPolyHandle = impoly(get(imHandle,'Parent'), []);
        else
            userData.imPolyHandle = impoly(get(imHandle,'Parent'),userData.ROI);
        end
        setPositionConstraintFcn(userData.imPolyHandle,fcn);
    end
end

set(handles.figure1, 'UserData', userData);
guidata(hObject,handles);
impoly_finish_cb(hObject,eventdata,handles)

function impoly_finish_cb(hObject,eventdata,handles)
% impoly_finish_cb displays a useful help message if the user double clicks
% on the polygon selection
    userData = get(handles.figure1, 'UserData');
        while(isempty(userData.ROI) && userData.imPolyHandle.isvalid)
            resume(userData.imPolyHandle);
            wait(userData.imPolyHandle);
            % After user double clicks, check if figure is still valid
            if(isvalid(handles.figure1))
                userData = get(handles.figure1, 'UserData');
                % Check if a valid ROI has been set already, meaning we are closing
                % or saving. If not, then display a help message
                if(isempty(userData.ROI))
                    helpdlg(['Click "Save" to store the ROI, or ' ...
                             'Click "Draw new region of interest" to reset. ' ...
                             'See impoly documentation for other commands.']);
                end
            end
        end

function close_previewFig(hObject, eventdata)
handles = guidata(get(hObject,'UserData'));
userData=get(handles.figure1,'UserData');
if(userData.imPolyHandle.isvalid)
    userData.ROI=getPosition(userData.imPolyHandle);
    set(handles.figure1,'UserData',userData);
    update_data(hObject, eventdata, handles);
end

% --- Executes on slider movement.
function frameNumberEdition_Callback(hObject, eventdata, handles)
userData = get(handles.figure1, 'UserData');

% Retrieve the value of the selected image
if strcmp(get(hObject,'Tag'),'edit_frameNumber')
    frameNumber = str2double(get(handles.edit_frameNumber, 'String'));
else
    frameNumber = get(handles.slider_frameNumber, 'Value');
end
frameNumber=round(frameNumber);

% Check the validity of the frame values
if isnan(frameNumber)
    warndlg('Please provide a valid frame value.','Setting Error','modal');
end
frameNumber = min(max(frameNumber,1),userData.nFrames);

% Store value
set(handles.slider_frameNumber,'Value',frameNumber);
set(handles.edit_frameNumber,'String',frameNumber);

% Save data and update graphics
set(handles.figure1, 'UserData', userData);
guidata(hObject, handles);
update_data(hObject,eventdata,handles);

% --- Executes on button press in pushbutton_save.
function pushbutton_save_Callback(hObject, eventdata, handles)

userData = get(handles.figure1, 'UserData');

if get(handles.checkbox_addROI, 'Value')
    % Check valid output directory for the subregional movie
    outputDirectory = get(handles.edit_outputDirectory,'String');
    if isempty(outputDirectory),
        errordlg('Please select an output directory','Error','modal');
        return;
    end
end

% Read ROI if crop window is still visible
if userData.imPolyHandle.isvalid
    userData.ROI=getPosition(userData.imPolyHandle);
    set(handles.figure1,'UserData',userData);
end
update_data(hObject,eventdata,handles);

% Create ROI mask and save it in the outputDirectory
userData = get(handles.figure1, 'UserData');
try
    mask=createMask(userData.imPolyHandle);           
    
    if get(handles.checkbox_addROI, 'Value')
        
        %Save  mask to ROI output directory
        maskPath = fullfile(outputDirectory,'roiMask.tif');
        imwrite(mask,maskPath);
        
        % Create a new region of interest and save the object
        userData.MD.addROI(maskPath, outputDirectory);
        movieROI = userData.MD.rois_(end);
        movieROI.save();

        % If called from movieSelectorGUI
        if userData.mainFig ~=-1,
            % Retrieve main window userData
            userData_main = get(userData.mainFig, 'UserData');

            % Append new ROI to movie selector panel
            userData_main.MD = horzcat(userData_main.MD, movieROI);
            set(userData.mainFig, 'UserData', userData_main)
            movieSelectorGUI('refreshDisplay',userData.mainFig,...
                eventdata,guidata(userData.mainFig));
        end
    else
        %Save mask to this MD's output directory
        maskPath = fullfile(userData.MD.outputDirectory_,'roiMask.tif');
        imwrite(mask,maskPath);
        
        % Create a new region of interest and save the object
        userData.MD.setROIMaskPath(maskPath);
        userData.MD.save();
    end
    % Delete current window
    delete(handles.figure1)
catch err
    errordlg('Did not successfully save ROI. Complete ROI selection or Cancel.');
    rethrow(err);
end


% --- Executes on button press in pushbutton_outputDirectory.
function pushbutton_outputDirectory_Callback(hObject, eventdata, handles)

userData = get(handles.figure1, 'UserData');
pathname = uigetdir(userData.MD.getPath(),'Select output directory');

% Test uigetdir output and store its results
if isequal(pathname,0), return; end
set(handles.edit_outputDirectory,'String',pathname);

% Save data
set(handles.figure1,'UserData',userData);
guidata(hObject, handles);

% --- Executes on button press in checkbox_addROI.
function checkbox_addROI_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_addROI (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if get(hObject,'Value')
    state = 'on';
else
    state = 'off';
end
set(handles.pushbutton_outputDirectory, 'Enable', state);
set(handles.edit_outputDirectory, 'Enable', state);


% --- Executes on button press in pushbutton_draw.
function pushbutton_draw_Callback(hObject, eventdata, handles)

userData=get(handles.figure1,'UserData');
userData.ROI = [];
if userData.imPolyHandle.isvalid
    delete(userData.imPolyHandle);
end
set(hObject, 'UserData', userData);
update_data(hObject,eventdata,handles);
