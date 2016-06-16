function varargout = wobble_correct(varargin)
% WOBBLE_CORRECT MATLAB code for wobble_correct.fig
%      WOBBLE_CORRECT, by itself, creates a new WOBBLE_CORRECT or raises the existing
%      singleton*.
%
%      H = WOBBLE_CORRECT returns the handle to a new WOBBLE_CORRECT or the handle to
%      the existing singleton*.
%
%      WOBBLE_CORRECT('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in WOBBLE_CORRECT.M with the given input arguments.
%
%      WOBBLE_CORRECT('Property','Value',...) creates a new WOBBLE_CORRECT or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before wobble_correct_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to wobble_correct_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help wobble_correct

% Last Modified by GUIDE v2.5 12-Jun-2016 12:29:36

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @wobble_correct_OpeningFcn, ...
                   'gui_OutputFcn',  @wobble_correct_OutputFcn, ...
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


% --- Executes just before wobble_correct is made visible.
function wobble_correct_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to wobble_correct (see VARARGIN)

% Choose default command line output for wobble_correct
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes wobble_correct wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = wobble_correct_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in pushbutton1.
function pushbutton1_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% select localization file
[fname, fpath] = uigetfile('*.*');

fullname = fullfile(fpath,fname);

set(handles.text5,'String',fullname);


% --- Executes on button press in pushbutton2.
function pushbutton2_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% select ground truth file
[fname, fpath] = uigetfile('*.csv');
fullname = fullfile(fpath,fname);
set(handles.text6,'String', fullname);



function edit1_Callback(hObject, eventdata, handles)
% hObject    handle to edit1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit1 as text
%        str2double(get(hObject,'String')) returns contents of edit1 as a double
%if ~exist(get(hObject,'String'),'dir')
%    mkdir(get(hObject,'String'));
%end

% --- Executes during object creation, after setting all properties.
function edit1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton3.
function pushbutton3_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
foldpath = uigetdir(pwd,'Select a folder to store the output file');
 
set(handles.edit1,'String', foldpath);


% --- Executes on button press in pushbutton4.
function pushbutton4_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
fullnameLoc = get(handles.text5,'String');
%hasHeader = fgetl(fopen(fullnameLoc));
%hasHeader = 1*(sum(isstrprop(hasHeader,'digit'))/length(hasHeader) < .6);
%localData = csvread(fullnameLoc, hasHeader, 0);
%SH: switched to importdata tool and defined columns to make more general
localData =importdata(fullnameLoc);
if isstruct(localData)
    %strip the header
    localData = localData.data;
end
xCol = str2num(get(handles.edit_x,'String'));
yCol = str2num(get(handles.editY,'String'));
frCol = str2num(get(handles.editFr,'String'));

fullnameGT = get(handles.text6,'String');
%assumes GT file is as defined in competition
%CSV file. X col 3, y col 4.
gtData = importdata(fullnameGT);
XCOLGT =3;
YCOLGT =4;
gtAll = gtData(:,[XCOLGT,YCOLGT]);
gt = unique(gtAll,'rows');

frameIsOneIndexed = get(handles.radiobutton_is1indexed,'Value');

[pathstr,~,~] = fileparts(fullnameLoc); 
output_path = pathstr;
xnm = localData(:,xCol);
ynm = localData(:,yCol);
frame = localData(:,frCol);

%might be set by the users in future updates
zmin = -750;zmax = 750;zstep = 10;%nm
roiRadius = 500;%nm

wobbleCorrectSimBead(xnm,ynm,frame, gt,zmin,zstep,zmax,roiRadius,frameIsOneIndexed,output_path)


% --- Executes on button press in pushbutton5.
function pushbutton5_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)



function edit_x_Callback(hObject, eventdata, handles)
% hObject    handle to edit_x (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_x as text
%        str2double(get(hObject,'String')) returns contents of edit_x as a double


% --- Executes during object creation, after setting all properties.
function edit_x_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_x (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function editY_Callback(hObject, eventdata, handles)
% hObject    handle to editY (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editY as text
%        str2double(get(hObject,'String')) returns contents of editY as a double


% --- Executes during object creation, after setting all properties.
function editY_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editY (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function editFr_Callback(hObject, eventdata, handles)
% hObject    handle to editFr (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editFr as text
%        str2double(get(hObject,'String')) returns contents of editFr as a double


% --- Executes during object creation, after setting all properties.
function editFr_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editFr (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
%-------------------------------------------------------
function wobbleCorrectSimBead(xnm,ynm,frame,gt,zmin,zstep,zmax,roiRadius,frameIsOneIndexed,output_path)

nBead = size(gt,1);
for ii = 1:nBead
    beadPos = gt(ii,:);
    ROInm(ii,1:2) = beadPos-roiRadius;%xmin ymin
    ROInm(ii,3:4) = 2*roiRadius;%width height
end

zSlice = zmin:zstep:zmax;

%frameIsOneIndexed = ~sum(frame==0) > 0;%should detect it automatically
if ~frameIsOneIndexed
    frame = frame+1;%have to account for possilbe zero-indexing or everthing will get screwed up
end

znm = zSlice(frame)';

wobbleMatrix = wobbleCalibration(xnm, ynm, znm, nBead, 'ROI', ROInm, 'Zfit', znm, 'NumSplineBreak', 10,...
    'GT', gt);
[~, indCorr] = unique(wobbleMatrix(:,1));
wobbleMatrixUnique = wobbleMatrix(indCorr,[2,3,1]);
%save in csv file, units : nm, column order : X Y Z
csvwrite(fullfile(output_path,'wobbleCorrectionData.csv'), wobbleMatrixUnique);
saveas(gcf,fullfile(output_path,'XY wobble result.fig'));
saveas(gcf,fullfile(output_path,'XY wobble result.png'));

%------------------------------------------------------------------------

function [wobbleMatrix] = wobbleCalibration(x,y,z,nBead,varargin)
% WOBBLECALIBRATION Generate correction data for z-dependent "wobble"
%
% MODIFIED BY THANH-AN 24th May 2016
%
% SYNTAX
%  [wobbleMatrix]=wobbleCalibration(x,y,z,nBead)
%  [wobbleMatrix]=wobbleCalibration(..., 'SaveWobbleFile',wobbleFileName)
%  [wobbleMatrix]=wobbleCalibration(..., 'ZLimit',zLim)
%  [wobbleMatrix]=wobbleCalibration(..., 'NumSplineBreak',nBreak)
%
% INPUTS
%  x,y,z: A dataset containing the observed localizations 'x,y' of beads at known position 'z'
%     Usually generated by imaging fluorescent beads, stepped through z with a 
%     z-piezo stage. Normally, the exact same data used to generate an astigmatic
%     3D PSF width calibration is used.
%  nBead: The number of 'good' (ie non-aggregated beads) in the image, which will
%     be manually selected.
% OUTPUTS
%  wobbleMatrix
%     Lookup table containing a list of axial shifts 
%     as a function of Z, ie 
%     wobbleData = [Z1, xShift1, yShift1;...
%                   Z2, xShift2, yShift2;... etc]
% DESCRIPTION
%  [wobbleMatrix]=wobbleCalibration(x,y,z,nBead)
%     Plots the bead localizations and prompts the user to manually select 'nBead'
%     regions within the image containing only localizations from a single bead.
%     The xy-wobble is calculated as a function of z for each bead via spline fitting, and plotted.
%     The combined xy-wobble is calculated via a second spline fit to exclude outliers usually arising
%     from aggregated or overlapping beads.
%     A wobble lookup table 'wobbleMatrix' is generated, which may be used for subsequent correction
%     of 3D localization data.
%
% OPTIONS
%  [wobbleMatrix]=wobbleCalibration(..., 'SaveWobbleFile',wobbleFileName)
%     Save the 'wobbleMatrix' to a text file 'wobbleFileName'
%  [wobbleMatrix]=wobbleCalibration(..., 'ZLimit',zLim)
%     Only produce a wobble lookup table over specified limits 
%     zLim = [zMin zMax]. This is useful for excluding regions where fitting is unreliable 
%     because the beads have become too defocussed. Ie, bead calibration data is usually
%     taken over a large Z range. This allows cropping to only the useful z-range
%  [wobbleMatrix]=wobbleCalibration(..., 'NumSplineBreak',nBreak)
%     Set the number of breaks 'nBreak' (default = 10) to use in the spline fit. If the resolution of the 
%     spline is insufficient compared to the underlying data, consider increasing 'nBreak'
%  [wobbleMatrix]=wobbleCalibration(..., 'ROI',ROI)
%     Provides the beads ROIs instead of the manual selection
%  [wobbleMatrix]=wobbleCalibration(..., 'Zfit',zfit)
%     Set the z-slice for which the XY correction are calculated (e.g.
%     -750 nm to 750 nm every 10 nm.
%  [wobbleMatrix]=wobbleCalibration(..., 'GT',gt)
%     Calculate the xy-shift wrt the ground truth gt instead of the xy @z=0
%     It doesn't matter if the beads positions are unordered wrt ROI

%  EXAMPLE
%  Generate a wobble correction lookup table for the test data supplied with these functions.
%
%  The test data below ('bead 0.1um -1.5to1.5 20nm Z-step 3D cal.txt') shows fluorescent beads 
%  (Invitrogen 0.1um TetraSpeks), stepped in Z with a piezo stage and 2D localized with RapidSTORM 
%  in X,Y. Z was set to the known position of the piezo via RapidSTORM. This dataset can similarly
%  be used to generate a PSF width 3D lookup table for astigmatic Z-localization (see the 
%  RapidSTORM 3.3 manual for further instructions on how to do this). 
%
%  Load the x,y,z data (this is the test data supplied with the 
%  wobble correction functions):
%     fname = 'test data\bead 0.1um -1.5to1.5 20nm Z-step 3D cal.txt'
%     a=importdata(fname);data=a.data;
%     x = data(:,1); y=data(:,3);z=data(:,5);
%  You also need to tell the progam how many good beads are in the image. Do this by loading
%  up the localizations in your favorite PALM visualization software (eg PALMsiever
%  https://github.com/PALMsiever/palm-siever), and counting how many good beads you have
%     nBead = 7;
%  Run the calibration, saving the output
%     [wobbleMatrix]=wobbleCalibration(x,y,z,nBead,'SaveWobbleFile','Wobble-cal test.txt');
%  A scatter plot of the XY localizations will appear, you will be prompted to select 'nBead' 
%  (here, 7) rectangular bead-containing regions.
%  Once selected, a plot of XY-wobble vs z for each bead should be generated, together
%  with combined fits for all beads.
%  Note that the fit, and hence wobble correction, becomes unreliable once the beads go 
%  out of focus (here z<-750 and z>850). In practice, Z-localization in these regions is also 
%  unlikely to be feasible. Therefore, exclude these regions from the lookup table, 
%  either by manually editing the wobble file, or by re-running the calibration with:
%     [wobbleMatrix]=wobbleCalibration(x,y,z,nBead, ...,
%                       'SaveWobbleFile','Wobble-cal test.txt','ZLimit',[-750 850]);
%  The advantage of rerunning the calibration like this is that the (default) 10 spline points 
%  are spread over a smaller range, giving higher resolution to the spline fit. Alternatively
%  run the entire range with a higher number of spline points to begin with:
%     [wobbleMatrix]=wobbleCalibration(x,y,z,nBead, ...,
%                       'SaveWobbleFile','Wobble-cal test.txt','NumSplineBreak',20);
%  and manually crop the text file later.
%
%  The generated wobbleMatrix may now be used for wobble correction. See CORRECTWOBBLE documentation
%  for details.
%
% This software is released under the GPL v3 (see license file 'gpl.txt'). It is provided AS-IS and no
% warranty is given.
%
% Author: Seamus Holden
% Last update: April 2015


narg = numel(varargin);
nBreak = 10;
zLim = [-Inf Inf];
ii=1;
%doSaveFile = false;
hasROI = false;
wobbleSaveName = [];
zfit = [];
gt = [];

while ii<=narg
   if strcmp(varargin{ii},'NumSplineBreak')
      nBreak= varargin{ii+1};
      ii = ii+2;
   elseif strcmp(varargin{ii},'ZLimit')
      zLim= varargin{ii+1};
      ii = ii+2;
   elseif strcmp(varargin{ii},'SaveWobbleFile')
      wobbleSaveName= varargin{ii+1};
      ii = ii+2;
   elseif strcmp(varargin{ii},'ROI')
       hasROI = true;
       ROI = varargin{ii+1};
       ii = ii+2;
   elseif strcmp(varargin{ii},'Zfit')
       zfit = varargin{ii+1};
       ii = ii + 2;
   elseif strcmp(varargin{ii},'GT')
       gt = varargin{ii+1};
       ii = ii + 2;
   else
      ii = ii+1;
   end
end
%Modified by Thanh-an Pham the 16th May 2016


if hasROI
    for ii = 1:nBead
        bead{ii} = [ROI(ii,1),ROI(ii,2),ROI(ii,3) + ROI(ii,1), ROI(ii,4) + ROI(ii,2)];
    end
else
    figure;
    hF = scatter(x,y,25,z,'.');
    set(gca,'YDir','reverse')
    
    for ii = 1:nBead
        hR{ii} = imrect(hF);
        bead{ii} = getPosition(hR{ii});
        bead{ii}(3) = bead{ii}(3) + bead{ii}(1);
        bead{ii}(4) = bead{ii}(4) + bead{ii}(2);
    end
end


% bead lim are [xmin, ymin, xmax, ymax]

[z,xWobble, yWobble] = xyWobble(x,y,z,bead,zLim,wobbleSaveName,nBreak,zfit,gt);

wobbleMatrix = [z(:),xWobble(:),yWobble(:)];
%---------------------------------------------------
function [zfit,xWobble, yWobble] = xyWobble(x,y,z,beadLim,zlim,fsavename,nBreak,zfit,gt)
WOBBLEWARNINGNM = 300;%warn if values greater than this


bead = beadLim;

n = numel(bead);
zAll = [];
gt_tmp = gt;
k=1;
for ii = 1:n
   isBead = x>bead{ii}(1) & y>bead{ii}(2) & x<bead{ii}(3) & y<bead{ii}(4);
   xBead{ii} = x(isBead);
   yBead{ii} = y(isBead);
   zBead{ii} = z(isBead);
   %reorder ground truth
   for jj = 1:n
       if ~isempty(gt_tmp) &&...
               gt_tmp(jj,1) > bead{ii}(1) && gt_tmp(jj,2) > bead{ii}(2) &&...
               gt_tmp(jj,1) < bead{ii}(3) && gt_tmp(jj,2) < bead{ii}(4)
           gt(ii,:) = gt_tmp(jj,:);
       end
   end
   zAll = [zAll;z(isBead)];
end

isOk = zAll>zlim(1)&zAll<zlim(2);
zRangeSet = zAll(isOk);

%Modified by Thanh-an Pham 16.05.2016
if isempty(zfit)
    zfit = min(zRangeSet): (max(zRangeSet)-min(zRangeSet))/nBreak:max(zRangeSet);
    %zfit = -750:10:750;
end
for ii =1:n
   xWobble = fit1Spline(xBead{ii},zBead{ii},zfit,nBreak);
   yWobble = fit1Spline(yBead{ii},zBead{ii},zfit,nBreak);
   beadFit{ii} = [zfit(:), xWobble(:),yWobble(:)];
end

if isempty(gt)
    %find the zfit point nearest to zero, align everything on this
    [~, idx] =min(abs(zfit));
    for ii =1:n
        %shift x
        beadFit{ii}(:,2) =  beadFit{ii}(:,2) - beadFit{ii}(idx,2);
        %shift y
        beadFit{ii}(:,3) =  beadFit{ii}(:,3) - beadFit{ii}(idx,3);
    end
else
    %use the ground truth gt for xy for each z
    for ii = 1:n
        %shift x
        beadFit{ii}(:,2) = beadFit{ii}(:,2) - gt(ii,1);
        %shift y
        beadFit{ii}(:,3) = beadFit{ii}(:,3) - gt(ii,2);
    end
end

%combine all the spline fits, one more spline fit to generate the final data
z=[];x=[];y=[];
for ii =1:n
   z= [z;beadFit{ii}(:,1)];
   x= [x;beadFit{ii}(:,2)];
   y= [y;beadFit{ii}(:,3)];
end
xWobble = fit1Spline(x,z,zfit,nBreak);
yWobble = fit1Spline(y,z,zfit,nBreak);


%plot
figure;hold all
plot(zfit,xWobble,'r');
plot(zfit,yWobble,'b');
for ii = 1:n
   plot(beadFit{ii}(:,1),beadFit{ii}(:,2),'k');
   plot(beadFit{ii}(:,1),beadFit{ii}(:,3),'g');
end
legend('X, combined fit','Y, combined fit', 'X, single bead fit', 'Y, single bead fit');
xlabel('Z (nm)');
ylabel('XY wobble (nm)')
%saveas(gcf,'XY wobble result.fig');
%saveas(gcf,'XY wobble result.png');

calData = [zfit(:), xWobble(:),yWobble(:)];
if ~isempty(fsavename)
   dlmwrite(fsavename, calData,' ');
end


if any([xWobble(:);yWobble]>=WOBBLEWARNINGNM)
    warning(['Wobble correction values > ', num2str(WOBBLEWARNINGNM),' nm detected, please check the input data for errors.']);
end
%-----------------------------------------
function [xfit] = fit1Spline(x,t,tfit,nBreak)

%fit with splinefit
ppX=splinefit(t,x,nBreak,'r');

xfit =  ppval(ppX,tfit);

%-----------------------------------------
function pp = splinefit(varargin)
%SPLINEFIT Fit a spline to noisy data.
%   PP = SPLINEFIT(X,Y,BREAKS) fits a piecewise cubic spline with breaks
%   (knots) BREAKS to the noisy data (X,Y). X is a vector and Y is a vector
%   or an ND array. If Y is an ND array, then X(j) and Y(:,...,:,j) are
%   matched. Use PPVAL to evaluate PP.
%
%   PP = SPLINEFIT(X,Y,P) where P is a positive integer interpolates the
%   breaks linearly from the sorted locations of X. P is the number of
%   spline pieces and P+1 is the number of breaks.
%
%   OPTIONAL INPUT
%   Argument places 4 to 8 are reserved for optional input.
%   These optional arguments can be given in any order:
%
%   PP = SPLINEFIT(...,'p') applies periodic boundary conditions to
%   the spline. The period length is MAX(BREAKS)-MIN(BREAKS).
%
%   PP = SPLINEFIT(...,'r') uses robust fitting to reduce the influence
%   from outlying data points. Three iterations of weighted least squares
%   are performed. Weights are computed from previous residuals.
%
%   PP = SPLINEFIT(...,BETA), where 0 < BETA < 1, sets the robust fitting
%   parameter BETA and activates robust fitting ('r' can be omitted).
%   Default is BETA = 1/2. BETA close to 0 gives all data equal weighting.
%   Increase BETA to reduce the influence from outlying data. BETA close
%   to 1 may cause instability or rank deficiency.
%
%   PP = SPLINEFIT(...,N) sets the spline order to N. Default is a cubic
%   spline with order N = 4. A spline with P pieces has P+N-1 degrees of
%   freedom. With periodic boundary conditions the degrees of freedom are
%   reduced to P.
%
%   PP = SPLINEFIT(...,CON) applies linear constraints to the spline.
%   CON is a structure with fields 'xc', 'yc' and 'cc':
%       'xc', x-locations (vector)
%       'yc', y-values (vector or ND array)
%       'cc', coefficients (matrix).
%
%   Constraints are linear combinations of derivatives of order 0 to N-2
%   according to
%
%     cc(1,j)*y(x) + cc(2,j)*y'(x) + ... = yc(:,...,:,j),  x = xc(j).
%
%   The maximum number of rows for 'cc' is N-1. If omitted or empty 'cc'
%   defaults to a single row of ones. Default for 'yc' is a zero array.
%
%   EXAMPLES
%
%       % Noisy data
%       x = linspace(0,2*pi,100);
%       y = sin(x) + 0.1*randn(size(x));
%       % Breaks
%       breaks = [0:5,2*pi];
%
%       % Fit a spline of order 5
%       pp = splinefit(x,y,breaks,5);
%
%       % Fit a spline of order 3 with periodic boundary conditions
%       pp = splinefit(x,y,breaks,3,'p');
%
%       % Constraints: y(0) = 0, y'(0) = 1 and y(3) + y"(3) = 0
%       xc = [0 0 3];
%       yc = [0 1 0];
%       cc = [1 0 1; 0 1 0; 0 0 1];
%       con = struct('xc',xc,'yc',yc,'cc',cc);
%
%       % Fit a cubic spline with 8 pieces and constraints
%       pp = splinefit(x,y,8,con);
%
%       % Fit a spline of order 6 with constraints and periodicity
%       pp = splinefit(x,y,breaks,con,6,'p');
%
%   See also SPLINE, PPVAL, PPDIFF, PPINT

%   Author: Jonas Lundgren <splinefit@gmail.com> 2010

%   2009-05-06  Original SPLINEFIT.
%   2010-06-23  New version of SPLINEFIT based on B-splines.
%   2010-09-01  Robust fitting scheme added.
%   2010-09-01  Support for data containing NaNs.
%   2011-07-01  Robust fitting parameter added.

% Check number of arguments
error(nargchk(3,7,nargin));

% Check arguments
[x,y,dim,breaks,n,periodic,beta,constr] = arguments(varargin{:});

% Evaluate B-splines
base = splinebase(breaks,n);
pieces = base.pieces;
A = ppval(base,x);

% Bin data
[junk,ibin] = histc(x,[-inf,breaks(2:end-1),inf]); %#ok

% Sparse system matrix
mx = numel(x);
ii = [ibin; ones(n-1,mx)];
ii = cumsum(ii,1);
jj = repmat(1:mx,n,1);
if periodic
    ii = mod(ii-1,pieces) + 1;
    A = sparse(ii,jj,A,pieces,mx);
else
    A = sparse(ii,jj,A,pieces+n-1,mx);
end

% Don't use the sparse solver for small problems
if pieces < 20*n/log(1.7*n)
    A = full(A);
end

% Solve
if isempty(constr)
    % Solve Min norm(u*A-y)
    u = lsqsolve(A,y,beta);
else
    % Evaluate constraints
    B = evalcon(base,constr,periodic);
    % Solve constraints
    [Z,u0] = solvecon(B,constr);
    % Solve Min norm(u*A-y), subject to u*B = yc
    y = y - u0*A;
    A = Z*A;
    v = lsqsolve(A,y,beta);
    u = u0 + v*Z;
end

% Periodic expansion of solution
if periodic
    jj = mod(0:pieces+n-2,pieces) + 1;
    u = u(:,jj);
end

% Compute polynomial coefficients
ii = [repmat(1:pieces,1,n); ones(n-1,n*pieces)];
ii = cumsum(ii,1);
jj = repmat(1:n*pieces,n,1);
C = sparse(ii,jj,base.coefs,pieces+n-1,n*pieces);
coefs = u*C;
coefs = reshape(coefs,[],n);

% Make piecewise polynomial
pp = mkpp(breaks,coefs,dim);


%--------------------------------------------------------------------------
function [x,y,dim,breaks,n,periodic,beta,constr] = arguments(varargin)
%ARGUMENTS Lengthy input checking
%   x           Noisy data x-locations (1 x mx)
%   y           Noisy data y-values (prod(dim) x mx)
%   dim         Leading dimensions of y
%   breaks      Breaks (1 x (pieces+1))
%   n           Spline order
%   periodic    True if periodic boundary conditions
%   beta        Robust fitting parameter, no robust fitting if beta = 0
%   constr      Constraint structure
%   constr.xc   x-locations (1 x nx)
%   constr.yc   y-values (prod(dim) x nx)
%   constr.cc   Coefficients (?? x nx)

% Reshape x-data
x = varargin{1};
mx = numel(x);
x = reshape(x,1,mx);

% Remove trailing singleton dimensions from y
y = varargin{2};
dim = size(y);
while numel(dim) > 1 && dim(end) == 1
    dim(end) = [];
end
my = dim(end);

% Leading dimensions of y
if numel(dim) > 1
    dim(end) = [];
else
    dim = 1;
end

% Reshape y-data
pdim = prod(dim);
y = reshape(y,pdim,my);

% Check data size
if mx ~= my
    mess = 'Last dimension of array y must equal length of vector x.';
    error('arguments:datasize',mess)
end

% Treat NaNs in x-data
inan = find(isnan(x));
if ~isempty(inan)
    x(inan) = [];
    y(:,inan) = [];
    mess = 'All data points with NaN as x-location will be ignored.';
    warning('arguments:nanx',mess)
end

% Treat NaNs in y-data
inan = find(any(isnan(y),1));
if ~isempty(inan)
    x(inan) = [];
    y(:,inan) = [];
    mess = 'All data points with NaN in their y-value will be ignored.';
    warning('arguments:nany',mess)
end

% Check number of data points
mx = numel(x);
if mx == 0
    error('arguments:nodata','There must be at least one data point.')
end

% Sort data
if any(diff(x) < 0)
    [x,isort] = sort(x);
    y = y(:,isort);
end

% Breaks
if isscalar(varargin{3})
    % Number of pieces
    p = varargin{3};
    if ~isreal(p) || ~isfinite(p) || p < 1 || fix(p) < p
        mess = 'Argument #3 must be a vector or a positive integer.';
        error('arguments:breaks1',mess)
    end
    if x(1) < x(end)
        % Interpolate breaks linearly from x-data
        dx = diff(x);
        ibreaks = linspace(1,mx,p+1);
        [junk,ibin] = histc(ibreaks,[0,2:mx-1,mx+1]); %#ok
        breaks = x(ibin) + dx(ibin).*(ibreaks-ibin);
    else
        breaks = x(1) + linspace(0,1,p+1);
    end
else
    % Vector of breaks
    breaks = reshape(varargin{3},1,[]);
    if isempty(breaks) || min(breaks) == max(breaks)
        mess = 'At least two unique breaks are required.';
        error('arguments:breaks2',mess);
    end
end

% Unique breaks
if any(diff(breaks) <= 0)
    breaks = unique(breaks);
end

% Optional input defaults
n = 4;                      % Cubic splines
periodic = false;           % No periodic boundaries
robust = false;             % No robust fitting scheme
beta = 0.5;                 % Robust fitting parameter
constr = [];                % No constraints

% Loop over optional arguments
for k = 4:nargin
    a = varargin{k};
    if ischar(a) && isscalar(a) && lower(a) == 'p'
        % Periodic conditions
        periodic = true;
    elseif ischar(a) && isscalar(a) && lower(a) == 'r'
        % Robust fitting scheme
        robust = true;
    elseif isreal(a) && isscalar(a) && isfinite(a) && a > 0 && a < 1
        % Robust fitting parameter
        beta = a;
        robust = true;
    elseif isreal(a) && isscalar(a) && isfinite(a) && a > 0 && fix(a) == a
        % Spline order
        n = a;
    elseif isstruct(a) && isscalar(a)
        % Constraint structure
        constr = a;
    else
        error('arguments:nonsense','Failed to interpret argument #%d.',k)
    end
end

% No robust fitting
if ~robust
    beta = 0;
end

% Check exterior data
h = diff(breaks);
xlim1 = breaks(1) - 0.01*h(1);
xlim2 = breaks(end) + 0.01*h(end);
if x(1) < xlim1 || x(end) > xlim2
    if periodic
        % Move data inside domain
        P = breaks(end) - breaks(1);
        x = mod(x-breaks(1),P) + breaks(1);
        % Sort
        [x,isort] = sort(x);
        y = y(:,isort);
    else
        mess = 'Some data points are outside the spline domain.';
        warning('arguments:exteriordata',mess)
    end
end

% Return
if isempty(constr)
    return
end

% Unpack constraints
xc = [];
yc = [];
cc = [];
names = fieldnames(constr);
for k = 1:numel(names)
    switch names{k}
        case {'xc'}
            xc = constr.xc;
        case {'yc'}
            yc = constr.yc;
        case {'cc'}
            cc = constr.cc;
        otherwise
            mess = 'Unknown field ''%s'' in constraint structure.';
            warning('arguments:unknownfield',mess,names{k})
    end
end

% Check xc
if isempty(xc)
    mess = 'Constraints contains no x-locations.';
    error('arguments:emptyxc',mess)
else
    nx = numel(xc);
    xc = reshape(xc,1,nx);
end

% Check yc
if isempty(yc)
    % Zero array
    yc = zeros(pdim,nx);
elseif numel(yc) == 1
    % Constant array
    yc = zeros(pdim,nx) + yc;
elseif numel(yc) ~= pdim*nx
    % Malformed array
    error('arguments:ycsize','Cannot reshape yc to size %dx%d.',pdim,nx)
else
    % Reshape array
    yc = reshape(yc,pdim,nx);
end

% Check cc
if isempty(cc)
    cc = ones(size(xc));
elseif numel(size(cc)) ~= 2
    error('arguments:ccsize1','Constraint coefficients cc must be 2D.')
elseif size(cc,2) ~= nx
    mess = 'Last dimension of cc must equal length of xc.';
    error('arguments:ccsize2',mess)
end

% Check high order derivatives
if size(cc,1) >= n
    if any(any(cc(n:end,:)))
        mess = 'Constraints involve derivatives of order %d or larger.';
        error('arguments:difforder',mess,n-1)
    end
    cc = cc(1:n-1,:);
end

% Check exterior constraints
if min(xc) < xlim1 || max(xc) > xlim2
    if periodic
        % Move constraints inside domain
        P = breaks(end) - breaks(1);
        xc = mod(xc-breaks(1),P) + breaks(1);
    else
        mess = 'Some constraints are outside the spline domain.';
        warning('arguments:exteriorconstr',mess)
    end
end

% Pack constraints
constr = struct('xc',xc,'yc',yc,'cc',cc);


%--------------------------------------------------------------------------
function pp = splinebase(breaks,n)
%SPLINEBASE Generate B-spline base PP of order N for breaks BREAKS

breaks = breaks(:);     % Breaks
breaks0 = breaks';      % Initial breaks
h = diff(breaks);       % Spacing
pieces = numel(h);      % Number of pieces
deg = n - 1;            % Polynomial degree

% Extend breaks periodically
if deg > 0
    if deg <= pieces
        hcopy = h;
    else
        hcopy = repmat(h,ceil(deg/pieces),1);
    end
    % to the left
    hl = hcopy(end:-1:end-deg+1);
    bl = breaks(1) - cumsum(hl);
    % and to the right
    hr = hcopy(1:deg);
    br = breaks(end) + cumsum(hr);
    % Add breaks
    breaks = [bl(deg:-1:1); breaks; br];
    h = diff(breaks);
    pieces = numel(h);
end

% Initiate polynomial coefficients
coefs = zeros(n*pieces,n);
coefs(1:n:end,1) = 1;

% Expand h
ii = [1:pieces; ones(deg,pieces)];
ii = cumsum(ii,1);
ii = min(ii,pieces);
H = h(ii(:));

% Recursive generation of B-splines
for k = 2:n
    % Antiderivatives of splines
    for j = 1:k-1
        coefs(:,j) = coefs(:,j).*H/(k-j);
    end
    Q = sum(coefs,2);
    Q = reshape(Q,n,pieces);
    Q = cumsum(Q,1);
    c0 = [zeros(1,pieces); Q(1:deg,:)];
    coefs(:,k) = c0(:);
    % Normalize antiderivatives by max value
    fmax = repmat(Q(n,:),n,1);
    fmax = fmax(:);
    for j = 1:k
        coefs(:,j) = coefs(:,j)./fmax;
    end
    % Diff of adjacent antiderivatives
    coefs(1:end-deg,1:k) = coefs(1:end-deg,1:k) - coefs(n:end,1:k);
    coefs(1:n:end,k) = 0;
end

% Scale coefficients
scale = ones(size(H));
for k = 1:n-1
    scale = scale./H;
    coefs(:,n-k) = scale.*coefs(:,n-k);
end

% Reduce number of pieces
pieces = pieces - 2*deg;

% Sort coefficients by interval number
ii = [n*(1:pieces); deg*ones(deg,pieces)];
ii = cumsum(ii,1);
coefs = coefs(ii(:),:);

% Make piecewise polynomial
pp = mkpp(breaks0,coefs,n);


%--------------------------------------------------------------------------
function B = evalcon(base,constr,periodic)
%EVALCON Evaluate linear constraints

% Unpack structures
breaks = base.breaks;
pieces = base.pieces;
n = base.order;
xc = constr.xc;
cc = constr.cc;

% Bin data
[junk,ibin] = histc(xc,[-inf,breaks(2:end-1),inf]); %#ok

% Evaluate constraints
nx = numel(xc);
B0 = zeros(n,nx);
for k = 1:size(cc,1)
    if any(cc(k,:))
        B0 = B0 + repmat(cc(k,:),n,1).*ppval(base,xc);
    end
    % Differentiate base
    coefs = base.coefs(:,1:n-k);
    for j = 1:n-k-1
        coefs(:,j) = (n-k-j+1)*coefs(:,j);
    end
    base.coefs = coefs;
    base.order = n-k;
end

% Sparse output
ii = [ibin; ones(n-1,nx)];
ii = cumsum(ii,1);
jj = repmat(1:nx,n,1);
if periodic
    ii = mod(ii-1,pieces) + 1;
    B = sparse(ii,jj,B0,pieces,nx);
else
    B = sparse(ii,jj,B0,pieces+n-1,nx);
end


%--------------------------------------------------------------------------
function [Z,u0] = solvecon(B,constr)
%SOLVECON Find a particular solution u0 and null space Z (Z*B = 0)
%         for constraint equation u*B = yc.

yc = constr.yc;
tol = 1000*eps;

% Remove blank rows
ii = any(B,2);
B2 = full(B(ii,:));

% Null space of B2
if isempty(B2)
    Z2 = [];
else
    % QR decomposition with column permutation
    [Q,R,dummy] = qr(B2); %#ok
    R = abs(R);
    jj = all(R < R(1)*tol, 2);
    Z2 = Q(:,jj)';
end

% Sizes
[m,ncon] = size(B);
m2 = size(B2,1);
nz = size(Z2,1);

% Sparse null space of B
Z = sparse(nz+1:nz+m-m2,find(~ii),1,nz+m-m2,m);
Z(1:nz,ii) = Z2;

% Warning rank deficient
if nz + ncon > m2
	mess = 'Rank deficient constraints, rank = %d.';
	warning('solvecon:deficient',mess,m2-nz);
end

% Particular solution
u0 = zeros(size(yc,1),m);
if any(yc(:))
    % Non-homogeneous case
	u0(:,ii) = yc/B2;
    % Check solution
	if norm(u0*B - yc,'fro') > norm(yc,'fro')*tol
        mess = 'Inconsistent constraints. No solution within tolerance.';
        error('solvecon:inconsistent',mess)
	end
end


%--------------------------------------------------------------------------
function u = lsqsolve(A,y,beta)
%LSQSOLVE Solve Min norm(u*A-y)

% Avoid sparse-complex limitations
if issparse(A) && ~isreal(y)
    A = full(A);
end

% Solution
u = y/A;

% Robust fitting
if beta > 0
    [m,n] = size(y);
    alpha = 0.5*beta/(1-beta)/m;
    for k = 1:3
        % Residual
        r = u*A - y;
        rr = r.*conj(r);
        rrmean = sum(rr,2)/n;
        rrmean(~rrmean) = 1;
        rrhat = (alpha./rrmean)'*rr;
        % Weights
        w = exp(-rrhat);
        spw = spdiags(w',0,n,n);
        % Solve weighted problem
        u = (y*spw)/(A*spw);
    end
end


%-----------------------------------------
function qq = ppdiff(pp,j)
%PPDIFF Differentiate piecewise polynomial.
%   QQ = PPDIFF(PP,J) returns the J:th derivative of a piecewise
%   polynomial PP. PP must be on the form evaluated by PPVAL. QQ is a
%   piecewise polynomial on the same form. Default value for J is 1.
%
%   Example:
%       x = linspace(-pi,pi,9);
%       y = sin(x);
%       pp = spline(x,y);
%       qq = ppdiff(pp);
%       xx = linspace(-pi,pi,201);
%       plot(xx,cos(xx),'b',xx,ppval(qq,xx),'r')
%
%   See also PPVAL, SPLINE, SPLINEFIT, PPINT

%   Author: Jonas Lundgren <splinefit@gmail.com> 2009

if nargin < 1, help ppdiff, return, end
if nargin < 2, j = 1; end

% Check diff order
if ~isreal(j) || mod(j,1) || j < 0
    msgid = 'PPDIFF:DiffOrder';
    message = 'Order of derivative must be a non-negative integer!';
    error(msgid,message)
end

% Get coefficients
coefs = pp.coefs;
[m n] = size(coefs);

if j == 0
    % Do nothing
elseif j < n
    % Derivative of order J
    D = [n-j:-1:1; ones(j-1,n-j)];
    D = cumsum(D,1);
    D = prod(D,1);
    coefs = coefs(:,1:n-j);
    for k = 1:n-j
        coefs(:,k) = D(k)*coefs(:,k);
    end
else
    % Derivative kills PP
    coefs = zeros(m,1);
end

% Set output
qq = pp;
qq.coefs = coefs;
qq.order = size(coefs,2);

%-----------------------------------------
function output = ppint(pp,a,b)
%PPINT Integrate piecewise polynomial.
%   QQ = PPINT(PP,A) returns the indefinite integral from A to X of a
%   piecewise polynomial PP. PP must be on the form evaluated by PPVAL.
%   QQ is a piecewise polynomial on the same form. Default value for A is
%   the leftmost break of PP.
%
%   I = PPINT(PP,A,B) returns the definite integral from A to B.
%
%   Example:
%       x = linspace(-pi,pi,7);
%       y = sin(x);
%       pp = spline(x,y);
%       I = ppint(pp,0,pi)
%
%       qq = ppint(pp,pi/2);
%       xx = linspace(-pi,pi,201);
%       plot(xx,-cos(xx),xx,ppval(qq,xx),'r')
%
%   See also PPVAL, SPLINE, SPLINEFIT, PPDIFF

%   Author: Jonas Lundgren <splinefit@gmail.com> 2009

if nargin < 1, help ppint, return, end
if nargin < 2, a = pp.breaks(1); end

% Get coefficients and breaks
coefs = pp.coefs;
[m n] = size(coefs);
xb = pp.breaks;
pdim = prod(pp.dim);

% Interval lengths
hb = diff(xb);
hb = repmat(hb,pdim,1);
hb = hb(:);

% Integration
coefs(:,1) = coefs(:,1)/n;
y = coefs(:,1).*hb;
for k = 2:n
    coefs(:,k) = coefs(:,k)/(n-k+1);
    y = (y + coefs(:,k)).*hb;
end
y = reshape(y,pdim,[]);
I = cumsum(y,2);
I = I(:);
coefs(:,n+1) = [zeros(pdim,1); I(1:m-pdim)];

% Set preliminary indefinite integral
qq = pp;
qq.coefs = coefs;
qq.order = n+1;

% Set output
if nargin < 3
    % Indefinite integral from a to x
    if a ~= xb(1)
        I0 = ppval(qq,a);
        I0 = I0(:);
        I0 = repmat(I0,m/pdim,1);
        qq.coefs(:,n+1) = qq.coefs(:,n+1) - I0;
    end
    output = qq;
else
    % Definite integral from a to b
    output = ppval(qq,b) - ppval(qq,a);
end


%-----------------------------------------
