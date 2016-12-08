classdef CompareToGroundTruth<interfaces.DialogProcessor
    methods
        function obj=CompareToGroundTruth(varargin)        
                obj@interfaces.DialogProcessor(varargin{:});
%             obj.inputParameters={'numberOfLayers','sr_layerson','mainfile','cam_pixelsize_nm'};
%             end   
        end
        
        function out=run(obj,p)
           
        end
        function pard=guidef(obj)
            pard=guidef(obj);
        end
    end
end


function wobble_callback(a,b,obj)
filename=(obj.locData.files.file(1).name);
path=fileparts(filename);
gtfile=[path filesep 'activations.csv','load ground truth for data'];
if ~exist(gtfile,'file')
    ind=strfind(path,filesep); 
    gtfile=[path(1:ind(end))  'activations.csv'];
end
 
if ~exist(gtfile,'file')
    [fi,pa]=uigetfile(gtfile);
    if fi
    gtfile=[pa fi];
    else
        return
    end
end

gtData=csvread(gtfile);
% fullnameLoc = get(handles.text5,'String');
%hasHeader = fgetl(fopen(fullnameLoc));
%hasHeader = 1*(sum(isstrprop(hasHeader,'digit'))/length(hasHeader) < .6);
%localData = csvread(fullnameLoc, hasHeader, 0);
%SH: switched to importdata tool and defined columns to make more general
% localData =importdata(fullnameLoc);
% if isstruct(localData)
%     %strip the header
%     localData = localData.data;
% end
% xCol = str2num(get(handles.edit_x,'String'));
% yCol = str2num(get(handles.editY,'String'));
% frCol = str2num(get(handles.editFr,'String'));
% 
% fullnameGT = get(handles.text6,'String');
% %assumes GT file is as defined in competition
% %CSV file. X col 3, y col 4.
% gtData = importdata(fullnameGT);
XCOLGT =3;
YCOLGT =4;
gtAll = gtData(:,[XCOLGT,YCOLGT]);
gt = unique(gtAll,'rows');

% frameIsOneIndexed = get(handles.radiobutton_is1indexed,'Value');
% 
% [pathstr,~,~] = fileparts(fullnameLoc); 
% output_path = pathstr;
% xnm = localData(:,xCol);
% ynm = localData(:,yCol);
% frame = localData(:,frCol);
cam_pixelsize_nm=obj.getPar('cam_pixelsize_nm');
p=obj.getGuiParameters;
if p.shiftpix
shiftx=-0.5*cam_pixelsize_nm;
shifty=-0.5*cam_pixelsize_nm;
else
    shiftx=0;
    shifty=0;
end
%might be set by the users in future updates
zmin = -750;zmax = 750;zstep = 10;%nm
roiRadius = 500;%nm
frameIsOneIndexed=true;
output_path=path;
wobbleCorrectSimBead(double(obj.locData.loc.xnm+shiftx),double(obj.locData.loc.ynm+shifty),double(obj.locData.loc.frame), gt,zmin,zstep,zmax,roiRadius,frameIsOneIndexed,filename)

% addpath('External/SMLMChallenge')
% wobble_correct;
end

function pard=guidef(obj)

pard.t1.object=struct('Style','text','String','Reference');
pard.t1.position=[1,1.5];
pard.t2.object=struct('Style','text','String','Target');
pard.t2.position=[1,2.5];

pard.t3.object=struct('Style','text','String','search radius (nm)');
pard.t3.position=[1,3.5];
pard.t3.Width=1;
pard.searchradius.object=struct('Style','edit','String','100');
pard.searchradius.position=[1,4.5];
pard.searchradius.Width=.5;
% pard.onlyfiltered.object=struct('Style','checkbox','String','Export filtered (displayed) localizations.','Value',1);
% pard.onlyfiltered.position=[2,1];
% pard.onlyfiltered.Width=2;
% pard.t1.Width=2;


pard.checklayer.object=struct('Style','checkbox','String','','Value',1);
pard.checklayer.position=[2,1];
pard.checklayer.Width=.5;

pard.reflayer.object=struct('Style','popupmenu','String',{{'layer1','layer2','layer3','layer4','layer5'}},'Value',1);
pard.reflayer.position=[2,1.5];
pard.reflayer.Width=1;

pard.targetlayer.object=struct('Style','popupmenu','String',{{'layer1','layer2','layer3','layer4','layer5'}},'Value',2);
pard.targetlayer.position=[2,2.5];
pard.targetlayer.Width=1;

pard.checkdata.object=struct('Style','checkbox','String','','Value',0);
pard.checkdata.position=[3,1];
pard.checkdata.Width=.5;

pard.refdata.object=struct('Style','popupmenu','String',{{'layer1','layer2','layer3','layer4','layer5'}},'Value',1);
pard.refdata.position=[3,1.5];
pard.refdata.Width=1;

pard.targetdata.object=struct('Style','popupmenu','String',{{'layer1','layer2','layer3','layer4','layer5'}},'Value',2);
pard.targetdata.position=[3,2.5];
pard.targetdata.Width=1;


% pard.forceungrouped.object=struct('Style','checkbox','String','Force ungrouped localiyations','Value',1);
% pard.forceungrouped.position=[3,3];
% pard.forceungrouped.Width=2;
% 
% pard.shiftpix.object=struct('Style','checkbox','String','Shift by 0.5 camera pixels','Value',1);
% pard.shiftpix.position=[3,1];
% pard.shiftpix.Width=2;
% 
% 
% 
% 
% pard.shiftframe.object=struct('Style','checkbox','String','Shift frame by +1','Value',1);
% pard.shiftframe.position=[4,1];
% pard.shiftframe.Width=4;
% 
% pard.comparer.object=struct('Style','popupmenu','String',{{'2D, 2013','3D, 2016'}},'Value',2);
% pard.comparer.position=[5,1];
% pard.comparer.Width=4;
% 
% pard.wobblebutton.object=struct('Style','pushbutton','String','Wobble.m','Callback',{{@wobble_callback,obj}});
% pard.wobblebutton.position=[4,4];
% pard.wobblebutton.Width=1;
pard.syncParameters={{'filelist_short','refdata',{'String'}},{'filelist_short','targetdata',{'String'}}};
pard.plugininfo.type='ProcessorPlugin';
end