classdef CameraSelector<handle
    %UNTITLED2 Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        handle
        imloader
        cameras
        guihandles
        currentcam=1;
        currentstate=1;
    end
    
    methods
        function obj=CameraSelector
            makeGui(obj)
        end
        function loadimages(obj,file)      
                obj.imloader=imageloaderAll(file);
                [par,cam]=getCameraCalibration(obj.imloader);
                if isempty(cam)
                    answ=questdlg('camera not found. Create new camera?');
                    if strcmp(answ,'Yes')
                        createnewcamera(obj);
                    end
                end
        end
    end
    
end

function makeGui(obj)
width=800;
height=500;
lineheight=25;
posbutton=width*.8;
buttonwidth=width*.15;
if isempty(obj.handle)||~isvalid(obj.handle)
     obj.handle=figure('Units','normalized','Units','pixels','Position',[150,200,width,height],'Name','CameraSelector','NumberTitle','off');
     delete(obj.handle.Children);
end
%load images
hp=uicontrol('Style','pushbutton','String','Load images','Position',[posbutton height-50,buttonwidth,lineheight],'Callback',{@loadimages,obj});

hp=uicontrol('Style','pushbutton','String','Save','Position',[posbutton height-80,buttonwidth,lineheight],'Callback',{@savecameras,obj});


tcam=uitable(obj.handle,'Position',[10 height-lineheight*3.5 width-200,lineheight*3]);
tcam.ColumnName={'Camera Name','ID field','ID'};
tcam.Data={'Cam1','Cam_ID','001'};
wh=tcam.Position(3);
tcam.ColumnWidth={'auto',wh*.5,wh*.3};
tcam.CellSelectionCallback={@cellselect_cam,obj,2,3};
tcam.ColumnEditable=[ true false true];

tpar=uitable(obj.handle,'Position',[10 height-lineheight*13-10 width-40,lineheight*9.5]);
tpar.ColumnName={'Parameter','mode','fixvalue','metafield','Value','conversion','Converted'};

% parnames={'EMon','pixelsize','conversion','emgain','offset','roi','exposure','timediff','comment'};
% dat=cell(length(parnames),7);
% dat(:,1)=parnames;
% dat{1,2}='fix';
dat=intpartable;


wh=tpar.Position(3);
tpar.ColumnWidth={'auto','auto','auto',wh*.25,'auto',wh*.2,'auto'};
tpar.CellSelectionCallback={@cellselect,obj,'par'};
tpar.CellEditCallback={@celledit,obj};
tpar.ColumnFormat={'char',{'fix','metadata','state dependent'},'char','char','char','char','char'};
tpar.ColumnEditable=[false true true false false true false];
tpar.Data=dat;


tstates=uicontrol('Style','listbox','String','State 1','Position',[10 30 width*.2,lineheight*4],'Callback',{@statecallback,obj});
tstatesadd=uicontrol('Style','pushbutton','String','add','Position',[10 lineheight*3+60 width*.08,lineheight],'Callback',{@stateadd,obj,1});
tstatesrem=uicontrol('Style','pushbutton','String','rem','Position',[width*.12 lineheight*3+60 width*0.08,lineheight],'Callback',{@stateadd,obj,2});
uicontrol('Style','text','String','State defining parameters','Position',[width*.25 lineheight*3+60 width*.2,lineheight])

tdef=uitable(obj.handle,'Position',[width*.25 30 width*.25,lineheight*4]);
tdef.ColumnName={'Meta Field','Value'};
tdef.Data={'select','';'select','';'select','';'select',''};
tdef.CellSelectionCallback={@cellselect,obj,'def'};

tval=uitable(obj.handle,'Position',[width*.7 30 width*.25,lineheight*5]);
tval.ColumnName={'Parameter','Value'};
tval.ColumnEditable=[false true];
tval.Data=tpar.Data(:,[1 3]);



obj.guihandles.camtable=tcam;
obj.guihandles.partable=tpar;
obj.guihandles.statevaltable=tval;


showpartable(obj);

obj.cameras(1).par=tpar.Data;
obj.cameras(1).ID=struct('name',tcam.Data{1,1},'tag',tcam.Data{1,2},'value',tcam.Data{1,3});
loadcameras(obj);
end

function celledit(table,data,obj)
end

function showpartable(obj)
partable=obj.guihandles.partable.Data;
tval=obj.guihandles.statevaltable.Data;

dat=(partable(:,[1 2]));

col=ones(size(tval,1),3)*.7;
indg=strcmp(dat(:,2),'state dependent');
col(indg,:)=1;
obj.guihandles.statevaltable.BackgroundColor=col;
obj.guihandles.statevaltable.Data=dat;
end

function loadcameras(obj)
file='settings/cameras.mat';
l=load(file);
obj.cameras=l.cameras;
obj.guihandles.camtable.Data=l.camtab;
obj.guihandles.partable.Data=l.cameras(1).par;
end

function savecameras(a,b,obj)
file='settings/cameras.mat';
cameras=obj.cameras;
cameras(obj.currentcam).par=obj.guihandles.partable.Data;
camtab=obj.guihandles.camtable.Data;
s=size(camtab);
for k=1:s(1)
    cameras(k).ID=struct('name',camtab{k,1},'tag',camtab{k,2},'value',camtab{k,3});
end
save(file,'cameras','camtab')
end

function t=intpartable
t=cell(12,7);
parnames={'EMon','pixelsize','conversion','emgain','offset','roi','exposure','timediff','comment','numberOfFrames','Width','Height'};
mode={'metadata','fix','state dependent','metadata','state dependent','metadata','metadata','metadata','metadata','metadata','metadata','metadata'};
default={'true','0.1','1','100','100','','1','1','settings not initialized','0','0','0'};
conversion={'str2double(X)','str2double(X)','str2double(X)','str2double(X)','str2double(X)','str2num(X)','str2double(X)','str2double(X)','','str2double(X)','str2double(X)','str2double(X)'};
metafield={'select','select','select','select','select','select','select','select','select','select','select','select','select'};

for k=1:size(t,1)
    t{k,1}=parnames{k};
    t{k,2}=mode{k};
    t{k,3}=default{k};
    t{k,4}=metafield{k};
    t{k,6}=conversion{k};
end
end

function cellselect_cam(table,data,obj,indtag,indval)
if isempty(data.Indices)
    return
end
if data.Indices(2)==indtag
    ma=obj.imloader.getmetadatatags;
     tag = gettag(ma);
     if ~isempty(tag)
        table.Data(data.Indices(1),[indtag indval])=tag;
     end
end
obj.cameras(obj.currentcam).par=obj.guihandles.partable.Data;
obj.guihandles.partable.Data=obj.cameras(data.Indices(1)).par;
obj.currentcam=data.Indices(1);
end

function cellselect(table,data,obj,tname)
switch tname
    case 'par'
        indtag=4;
        indval=5;
        tab=1;
    case 'def'
        indtag=1;
        indval=2;
        tab=2;
end

if isempty(data.Indices)
    return
end

if data.Indices(2)==indtag
    ma=obj.imloader.getmetadatatags;
     tag = gettag(ma);
     if ~isempty(tag)
        table.Data(data.Indices(1),[indtag indval])=tag;
        X=tag{2};
        if size(table.Data,2)>2
            if ~isempty(table.Data{data.Indices(1),6})&&~isempty(X)
                X=eval(table.Data{data.Indices(1),6});
            end
        table.Data{data.Indices(1),7}=X;
        end
     end
end

end

function tag = gettag(ma)
f=figure;
tc=uitable(f);
[~,ind]=sortrows(ma(:,1));


tc.Data=ma(ind,:);
tc.Position(3)=f.Position(3)-30;
tc.Position(2)=100;
tc.CellSelectionCallback=@cellselecth;
w=tc.Position(3);
tc.ColumnWidth={w*.75,.2*w};
uicontrol('Style','pushbutton','String','Ok','Position',[200 10 100 20],'Callback',@buttoncallback)
uicontrol('Style','pushbutton','String','Cancel','Position',[10 10 100 20],'Callback',@buttoncallback)
pos=[];
waitfor(f)

    function buttoncallback(a,b)
        if ~isempty(pos)&&strcmp(a.String,'Ok')
        tag=tc.Data(pos(1),:);
        else
            tag={};
        end
        close(f)
    end
    function cellselecth(t,d)
        pos=d.Indices;
    end
end


function loadimages(a,b,obj)

[file path]=uigetfile('*.*');
if file
    obj.loadimages([path file]);
% obj.imloader=imageloaderAll([path file]);
end
end

function createnewcamera(obj)
l=length(obj.cameras)+1;
dat=obj.guihandles.camtable.Data;
dat(l,:)={'new','select',''};
obj.guihandles.camtable.Data=dat;
obj.cameras(l)=obj.cameras(1);
obj.cameras(l).par=intpartable;
obj.cameras(l).ID=struct('name','new','tag','select','value','');
cellselect_cam(obj.guihandles.camtable,struct('Indices',[l,1]),obj,0,0);
end

function statecallback(a,b,obj)
end

function stateadd(a,b,obj,addrem)
end