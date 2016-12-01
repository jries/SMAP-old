classdef CameraSelector<handle
    %UNTITLED2 Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        handle
        imloader
        cameras
        guihandles
        currentcam=1;
    end
    
    methods
        function obj=CameraSelector
            makeGui(obj)
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
hp=uicontrol('Style','pushbutton','String','Load images','Position',[posbutton height-20,buttonwidth,lineheight],'Callback',{@loadimages,obj});

hp=uicontrol('Style','pushbutton','String','Save','Position',[posbutton height-80,buttonwidth,lineheight],'Callback',{@savecameras,obj});


tcam=uitable(obj.handle,'Position',[10 height-lineheight*3 width-200,lineheight*2]);
tcam.ColumnName={'Camera Name','ID field','ID'};
tcam.Data={'Cam1','Cam_ID','001'};
wh=tcam.Position(3);
tcam.ColumnWidth={'auto',wh*.4,wh*.3};
tcam.CellSelectionCallback={@cellselect_cam,obj,2,3};
tcam.ColumnEditable=[ true false false];

tpar=uitable(obj.handle,'Position',[10 height-lineheight*12 width-40,lineheight*8]);
tpar.ColumnName={'Parameter','mode','fixvalue','metafield','Value','conversion','Converted'};

% parnames={'EMon','pixelsize','conversion','emgain','offset','roi','exposure','timediff','comment'};
% dat=cell(length(parnames),7);
% dat(:,1)=parnames;
% dat{1,2}='fix';
dat=intpartable;


wh=tpar.Position(3);
tpar.ColumnWidth={'auto','auto','auto',wh*.25,'auto',wh*.2,'auto'};
tpar.CellSelectionCallback={@cellselect,obj,4,5};
tpar.CellEditCallback={@celledit,obj};
tpar.ColumnFormat={'char',{'fix','metadata','state dependent'},'char','char','char','char','char'};
tpar.ColumnEditable=[false true true false false true false];
tpar.Data=dat;

obj.guihandles.camtable=tcam;
obj.guihandles.partable=tpar;
obj.cameras(1).par=tpar.Data;
obj.cameras(1).ID=struct('name',tcam.Data{1,1},'tag',tcam.Data{1,2},'value',tcam.Data{1,3});
loadcameras(obj);
end

function celledit(table,data,obj)
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
parnames={'EMon','pixelsize','conversion','emgain','offset','roi','exposure','timediff','comment','numberOfFrames','Widht','Height'};
mode={'metadata','fix','state dependent','metadata','state dependent','metadata','metadata','metadata','metadata','metadata','metadata','metadata'};
default={'true','0.1','1','100','100','','1','1','settings not initialized','0','0','0'};
conversion={'str2double(X)','str2double(X)','str2double(X)','str2double(X)','str2double(X)','str2num(X)','str2double(X)','str2double(X)','','str2double(X)','str2double(X)','str2double(X)'};
metafield={'select','select','select','select','select','select','select','select','select','select','select','select','select'};

for k=1:size(t,1);
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

function cellselect(table,data,obj,indtag,indval)
if isempty(data.Indices)
    return
end
if data.Indices(2)==indtag
    ma=obj.imloader.getmetadatatags;
     tag = gettag(ma);
     if ~isempty(tag)
        table.Data(data.Indices(1),[indtag indval])=tag;
        X=tag{2};
        if ~isempty(table.Data{data.Indices(1),6})&&~isempty(X)
        X=eval(table.Data{data.Indices(1),6});
        end
        table.Data{data.Indices(1),7}=X;
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
obj.imloader=imageloaderAll([path file]);
end
end