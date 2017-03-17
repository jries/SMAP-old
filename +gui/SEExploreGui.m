classdef SEExploreGui<interfaces.SEProcessor
    properties
        hlines
        infostruct
        anglehandle
        sideviewax
    end
    methods
        function obj=SEExploreGui(varargin)   
            obj@interfaces.SEProcessor(varargin{:})
        end
        function makeGui(obj)
%             obj.handle=figure(37);
%             clf
            if ispc
                fontsize=14;
            else
                        fontsize=16;
            end
            set(obj.handle,'Position',[560,300,800,1000]);
             set(obj.handle,'MenuBar','none','Toolbar','none','SizeChangedFcn',{@sizechanged_callback,obj},'Name','ROIManager','NumberTitle','off')
             h.siteax=axes('Position',[.05,.6,.4,.32],'DataAspectRatio',[1 1 1],'NextPlot','replacechildren');
             h.cellax=axes('Position',[.55,.6,.4,.32],'DataAspectRatio',[1 1 1],'NextPlot','replacechildren');
             h.fileax=axes('Position',[.05,.2,.4,.32],'DataAspectRatio',[1 1 1],'NextPlot','replacechildren');%,'PlotBoxAspectRatio',[1 1 1],'DataAspectRatioMode','manual','PlotBoxAspectRatioMode','manual');
%              h.filelist=uicontrol(obj.handle,'Position',[.05,.05,.4,.1],'Style','listbox','String','X://','Units','normalized')
             h.filelist=uicontrol(obj.handle,'Position',[10,10,580+20,150],'Style','listbox','String','empty ','Units','normalized','FontSize',fontsize,'Callback',{@filelist_callback,obj});
             h.sitelist=uicontrol(obj.handle,'Position',[400-10,160,190+30,380],'Style','listbox',...
                 'String','empty ','Units','normalized','FontSize',fontsize,'max',100,'Callback',{@sitelist_callback,obj});
             h.celllist=uicontrol(obj.handle,'Position',[600+20,360,190-20,180],'Style','listbox',...
                 'String',' empty','Units','normalized','FontSize',fontsize,...
                 'Callback',{@celllist_callback,obj});
             h.redrawsite=uicontrol(obj.handle,'Position',[30,925,100,40],'Style','pushbutton','String','redraw','Units','normalized','FontSize',fontsize,'Callback',{@redrawsite_callback,obj});
             h.addsite=uicontrol(obj.handle,'Position',[290,925,60,40],'Style','pushbutton','String','Add','Units','normalized','FontSize',fontsize,'Callback',{@addsite,obj});
             h.removesite=uicontrol(obj.handle,'Position',[400-10,540,90,30],'Style','pushbutton','String','Remove','Units','normalized','FontSize',fontsize,'Callback',{@removesite_callback,obj});
             h.toggleuse=uicontrol(obj.handle,'Position',[510+20,540,70,30],'Style','pushbutton','String','Use','Units','normalized','FontSize',fontsize,'Callback',{@toggleuse_callback,obj});
             
             
             h.redrawcell=uicontrol(obj.handle,'Position',[430,925,100,40],'Style','pushbutton','String','redraw','Units','normalized','FontSize',fontsize,'Callback',{@redrawcell_callback,obj});
             h.addcell=uicontrol(obj.handle,'Position',[650,925,60,40],'Style','pushbutton','String','Add','Units','normalized','FontSize',fontsize,'Callback',{@addcell,obj});
             h.removecell=uicontrol(obj.handle,'Position',[600+20,540,90,30],'Style','pushbutton','String','Remove','Units','normalized','FontSize',fontsize,'Callback',{@removecell_callback,obj});
             
             h.redrawfile=uicontrol(obj.handle,'Position',[30,525,100,40],'Style','pushbutton','String','redraw','Units','normalized','FontSize',fontsize,'Callback',{@redrawfile_callback,obj});
             h.fileax.ButtonDownFcn={@fileaxclick,obj};
             h.cellax.ButtonDownFcn={@cellaxclick,obj};
             h.siteax.ButtonDownFcn={@siteaxclick,obj};
             
             h.info=uicontrol(obj.handle,'Position',[600+20,10,190-20,340],'Style','listbox','String','info','FontSize',fontsize*.85,'Max',10,'Units','normalized');
             
             h.anglebutton=uicontrol(obj.handle,'Position',[150,925,60,40],'Style','pushbutton','String','Angle','Units','normalized','FontSize',fontsize,'Callback',{@anglebutton_callback,obj});
             h.angle=uicontrol(obj.handle,'Position',[210,925,60,40],'Style','edit','String','0','Units','normalized','FontSize',fontsize,'Callback',{@angle_callback,obj});
              
             obj.guihandles=h;
             obj.updateFilelist;
             obj.hlines.rotationpos=[];
             obj.hlines.line1=[];
             obj.hlines.line2=[];
             obj.hlines.sidemarker.rotationpos=[];
             obj.hlines.sidemarker.line1=[];
             obj.hlines.sidemarker.line2=[];
             
             
             
             obj.handle.WindowKeyPressFcn={@keypress,obj,0};
             obj.addSynchronization('filelist_long',[],[],@obj.updateFilelist);
     
        end
 

        function updateSitelist(obj)
            redraw_sitelist(obj);
        end
        function updateCelllist(obj)
            redraw_celllist(obj);
        end
        function updateFilelist(obj,a,b)
            
            if ~isvalid(obj.handle)
                return
            end

            if ishandle(obj.guihandles.filelist)
                obj.guihandles.filelist.String={obj.SE.files.name};
                obj.guihandles.filelist.Value=1;
            end
            infofile=['settings' filesep 'infostruct.txt'];
            obj.infostruct=getinfostruct(infofile);
            redraw_sitelist(obj)
            redraw_celllist(obj)
        end
        
        function redrawall(obj,onlysites)
            global SMAP_stopnow
            if SMAP_stopnow
                disp('STOP button is activated. Function not executed')
            end
            if nargin<2
                onlysites=false;
            end
            if ~onlysites
        files=obj.SE.files;
        for k=1:length(files)
            obj.guihandles.filelist.Value=k;
            obj.status(['redrawall: file ' num2str(k) ' of ' num2str(length(files))]);
%             notify(obj.SE.locData,'status',recgui.statusEvent(['redrawall: file ' num2str(k) ' of ' num2str(length(files))]))
            drawnow
            filenumber=files(k).ID;
            files(k).image=[];
            obj.SE.plotfile(filenumber,obj.guihandles.fileax);
            if SMAP_stopnow
                break
            end

        end
            
        
        cells=obj.SE.cells;
        for k=1:length(cells)
            obj.guihandles.cellist.Value=k;   
            obj.status(['redrawall: cell ' num2str(k) ' of ' num2str(length(cells))])
            drawnow
            cells(k).image=[];
            obj.SE.plotcell(cells(k),obj.guihandles.cellax,obj.guihandles.fileax);
            

        end
        obj.SE.currentcell=cells(k);
        end
        sites=obj.SE.sites;
        indselected=obj.getSingleGuiParameter('sitelist').Value;
        if length(indselected)>1 %only redraw selected
            disp('redrawing only selected sites')
        else
            indselected=1:length(sites);
        end
        for k=indselected
            
            obj.guihandles.sitelist.Value=k;   
            obj.status(['redrawall: site ' num2str(k) ' of ' num2str(length(sites))])
            drawnow
            sites(k).image=[];
            obj.SE.plotsite(sites(k),obj.guihandles.siteax,obj.guihandles.cellax);
            obj.SE.processors.eval.evaluate(sites(k));
            sites(k).image.composite=[];
            sites(k).image.layers=[];
            sites(k).image.image=single(sites(k).image.image);
            if SMAP_stopnow
                break
            end
        end
        obj.SE.currentsite=sites(k);
        obj.status(['redrawall: completed'])
% 
%         obj.SE.plotsite(site,obj.guihandles.siteax,obj.guihandles.cellax);
        end
        function clearall(obj)
            obj.SE.clear;
            
            lDfiles=obj.locData.files.file;
            filenames={obj.locData.files.file(:).name};
            
            for k=1:length(filenames) 
                   f=filenames{k};
                    obj.SE.addFile(f,k,lDfiles(k).info);
            end   
            
            obj.updateFilelist
            redraw_celllist(obj)
            redraw_sitelist(obj)
        end
        function lineannotation(obj,linenumber,hline)
            if nargin<3||isempty(hline)           
                obj.anglehandle{linenumber}=obj.getPar(['ROI_lineannotation_handle_' num2str(linenumber)]);
                hline=obj.anglehandle{linenumber};
            end
            plotline(obj,['line',num2str(linenumber)],hline);
        end
        function nextsite(obj,direction)
            v=obj.guihandles.sitelist.Value+direction;
            s=obj.guihandles.sitelist.String;
            vnew=max(1,min(v,length(s)));
            obj.SE.currentsite=obj.SE.sites(vnew);
            obj.guihandles.sitelist.Value=vnew;
            plotsite(obj,obj.SE.currentsite)
        end
    end
end


function sizechanged_callback(object, event, obj)
% uiwait(obj.handle,1)  
% f=object.Position(3)/obj.guiPar.width;
% if f~=1
% obj.resize(f);
% obj.guiPar.width=object.Position(3);
% end
end

function toggleuse_callback(data,action,obj)
 selected=obj.guihandles.sitelist.Value;
 newstate=~obj.SE.sites(selected(1)).annotation.use;
if length(selected)<=1

site=obj.SE.currentsite;
site.annotation.use=newstate;
else
    for k=1:length(selected)
       site=obj.SE.sites(selected(k));
        site.annotation.use=newstate;
    end
end
obj.updateSitelist;
end

function fileaxclick(data,action,obj)
pos=action.IntersectionPoint*1000;
if action.Button==3
    %move
else
    
    currentcell=interfaces.SEsites;
    currentcell.pos=pos;
    currentcell.ID=0;

    currentcell.info.filenumber=obj.SE.currentfile.ID;
    obj.SE.currentcell=currentcell;
    obj.SE.plotcell(obj.SE.currentcell,obj.guihandles.cellax,obj.guihandles.fileax);
end

end

function cellaxclick(data,action,obj)
pos=action.IntersectionPoint*1000;
 pos(3)=0;
if action.Button==3
    obj.SE.currentcell.pos=pos;
    obj.SE.currentcell.image=[];
    obj.SE.updateCell;
    obj.SE.plotcell(obj.SE.currentcell,obj.guihandles.cellax,obj.guihandles.fileax);
else
    
    currentsite=interfaces.SEsites;
    currentsite.pos=pos;     
    currentsite.info.cell=obj.SE.currentcell.ID;
    currentsite.info.filenumber=obj.SE.currentfile.ID;
%     currentsite.sePar=obj.SE.sePar;
%     currentsite.annotation.rotationangle=0;
%     currentsite.annotation.rotationpos=zeros(2);
    obj.SE.currentsite=currentsite;
    plotsite(obj,obj.SE.currentsite);
end

end

function siteaxclick(data,action,obj)
pos=action.IntersectionPoint*1000;
if action.Button==3
    alphaimage=obj.SE.currentsite.image.angle;
    if alphaimage~=0
        posr=rotatepos(pos,obj.SE.currentsite.pos,-alphaimage);
        pos=posr;
    end
    
    obj.SE.currentsite.pos(1:2)=pos(1:2);
        obj.SE.updateSite;
%     end
%     obj.SE.plotsite(obj.SE.currentsite,obj.guihandles.siteax,obj.guihandles.cellax);
    obj.SE.currentsite.image=[];
    plotsite(obj,obj.SE.currentsite);
end
end

function addcell(data,action,obj)
obj.SE.currentcell.ID=obj.SE.addCell(obj.SE.currentcell);
plotfile(obj,obj.SE.currentfile.ID);
redraw_celllist(obj);
end

function addsite(data,action,obj)
if obj.SE.currentcell.ID==0
    button=questdlg('Add cell as well?');
    if strcmpi(button,'Yes')
        addcell(0,0,obj);
        obj.SE.currentsite.info.cell=obj.SE.currentcell.ID;
    else
        return
    end
end
obj.SE.currentsite.ID=obj.SE.addSite(obj.SE.currentsite);
% obj.guihandles.sitelist.Value=length(obj.guihandles.sitelist.String);
redraw_sitelist(obj);
% obj.guihandles.sitelist.Value=length(obj.guihandles.sitelist.String);
plotcell(obj,obj.SE.currentcell);
end


function redraw_celllist(obj)
if ~isvalid(obj.handle)
%     obj.delete;
    return
end
cells=obj.SE.cells;

if ~isempty(cells)&&~isempty(cells(1).ID)
    for k=1:length(cells)
        s{k}=['C' num2str(cells(k).ID,'%02.0f') 'F' num2str(cells(k).info.filenumber,'%02.0f')];
    end
    obj.guihandles.celllist.String=s;
    if isempty(obj.SE.currentcell)
        obj.SE.currentcell=obj.SE.cells(1);
    end
    val=obj.SE.indexFromID(obj.SE.cells,obj.SE.currentcell.ID);
    if isempty(val)
        val=obj.guihandles.celllist.Value;
    end
    val=min(max(1,val),length(s));
    obj.guihandles.celllist.Value=val;
else
    obj.guihandles.celllist.Value=1;
    obj.guihandles.celllist.String='empty';
end
end

function redraw_sitelist(obj)
if ~isvalid(obj.handle)
%     obj.delete;
    return
end
sites=obj.SE.sites;
if ~isempty(sites)&&~isempty(sites(1).ID)
for k=1:length(sites)
    usesite='';
    if isfield(sites(k).annotation,'use')
        if sites(k).annotation.use
            usesite='+';
        else
            usesite='-';
        end
    end
        
    
    list=[num2str(sites(k).annotation.list1.value) num2str(sites(k).annotation.list2.value)...
        num2str(sites(k).annotation.list3.value) num2str(sites(k).annotation.list4.value)];
    sitename=[num2str(sites(k).indList,'%2.0f') '.S' num2str(sites(k).ID,'%02.0f') 'C' num2str(sites(k).info.cell,'%02.0f')...
        'F' num2str(sites(k).info.filenumber,'%02.0f') 'L' list usesite];
    sites(k).name=sitename;
    s{k}=sitename;
end

obj.guihandles.sitelist.String=s;
obj.guihandles.sitelist.Value=max(1,min(obj.guihandles.sitelist.Value,length(s)));

else
    obj.guihandles.sitelist.Value=1;
    obj.guihandles.sitelist.String='empty';
end
end

function celllist_callback(data,action,obj)
ind=data.Value;
obj.SE.currentcell=obj.SE.cells(ind);
plotcell(obj,obj.SE.currentcell);
end

function plotcell(obj,cell)
%set to first site of cell or create new site centered in cell
 vold=obj.SE.files(obj.guihandles.filelist.Value).ID;
 if vold~=cell.info.filenumber
%      file=obj.SE.files(cell.info.filenumber);
     
     plotfile(obj,cell.info.filenumber);
     obj.SE.currentfile=obj.SE.files(cell.info.filenumber);
 end

newcellind=obj.SE.indexFromID(obj.SE.cells,cell.ID);
if newcellind>0
obj.guihandles.celllist.Value=newcellind;
end

obj.SE.plotcell(cell,obj.guihandles.cellax,obj.guihandles.fileax);
end


function sitelist_callback(data,action,obj)
if length(data.Value)==1
ind=data.Value(1);
obj.SE.currentsite=obj.SE.sites(ind);
plotsite(obj,obj.SE.currentsite);

    
end
end

function plotsite(obj, site)
% site.sePar=obj.SE.sePar;
vold=min(length(obj.SE.cells),obj.guihandles.celllist.Value);
% cellid=obj.SE.cells(vold).ID;
if vold>0&&site.info.cell>0&&obj.SE.cells(vold).ID~=site.info.cell
    newcellind=obj.SE.indexFromID(obj.SE.cells,site.info.cell);
%     obj.guihandles.celllist.Value=newcellind;
    obj.SE.currentcell=obj.SE.cells(newcellind);
    plotcell(obj,obj.SE.currentcell);
%     obj.SE.plotcell(obj.SE.currentcell,obj.guihandles.cellax,obj.guihandles.fileax);
%     obj.SE.plotfile(obj.SE.currentcell.info.filenumber,obj.guihandles.fileax);
end
if obj.getPar('se_drawsideview')&&isfield(obj.locData.loc,'znm')
    if isempty(obj.sideviewax)||~isvalid(obj.sideviewax)
        f=figure;
        obj.sideviewax=gca;
        obj.sideviewax.NextPlot='replacechildren';
        obj.sideviewax.ButtonDownFcn={@sideview_click,obj};
    end
    axz=obj.sideviewax;
else
    axz=[];
end
        
obj.SE.plotsite(site,obj.guihandles.siteax,obj.guihandles.cellax,axz);
obj.guihandles.angle.String=num2str(pos2angle((obj.SE.currentsite.annotation.rotationpos.pos)),'%2.1f');
if sum(obj.SE.currentsite.annotation.rotationpos.pos(:).^2)~=0
    anglebutton_callback(0,0,obj)
end

% anotation line
if obj.SE.currentsite.annotation.line1.length>0
obj.lineannotation(1)
end
if obj.SE.currentsite.annotation.line2.length>0
obj.lineannotation(2)
end
%plot info
%update annotations
obj.SE.processors.annotation.sitechange(obj.SE.currentsite);
obj.SE.processors.eval.evaluate(obj.SE.currentsite);

if ~obj.getPar('se_keeptempimgs')
    obj.SE.currentsite.image.composite=[];
    obj.SE.currentsite.image.layers=[];
end
% l=obj.SE.currentsite.image.layers;
% for k=length(l):-1:1
%     obj.SE.currentsite.image.imax=l(k).images.finalImages.imax;
% end

% 
obj.SE.currentsite.image.image=single(obj.SE.currentsite.image.image);
%plot info
plotinfo(obj,site)
end

function sideview_click(hax,dat,obj)
site=obj.SE.currentsite;
pos=dat.IntersectionPoint;
site.pos(3)=pos(2)*1000;
site.image=[];
plotsite(obj,site)


end

function plotinfo(obj,site)
textwidth=20;
info=obj.infostruct;
s={'no info'};
for k=1:length(info)
    try
        val=eval(['site' info(k).field]);
        if length(val)>1
            fs='%.0f,';
        elseif isnumeric(val)&&abs(max(val))>100
            fs='%.0f';
        else
            
            fs=3;
        end
    catch ME
        
        val='n.d.';
%         disp('info file error or field not defined')
    end
    vals=num2str(val,fs);
    
    txt=info(k).text;
    if length(txt)<textwidth&&length(vals)<10;
        txt(end+1:textwidth)=' ';
    end
    
    s{k}=[txt vals];

end
obj.guihandles.info.String=s;
end



function filelist_callback(data,action,obj)
    filenumber=obj.SE.files(data.Value).ID;

obj.SE.currentfile=obj.SE.files(filenumber);
plotfile(obj,filenumber);
end

function plotfile(obj,filenumber)
obj.SE.plotfile(filenumber,obj.guihandles.fileax)
fileind=obj.SE.indexFromID(obj.SE.files,filenumber);
obj.guihandles.filelist.Value=fileind;
end


function redrawsite_callback(a,b,obj)
obj.SE.currentsite.image=[];
plotsite(obj,obj.SE.currentsite);
end

function redrawcell_callback(a,b,obj)
obj.SE.currentcell.image=[];
plotcell(obj,obj.SE.currentcell);
end
function redrawfile_callback(a,b,obj)
fileID=obj.SE.currentfile.ID;
indf=obj.SE.indexFromID(obj.SE.files,fileID);
obj.SE.files(indf).image=[];
plotfile(obj,obj.SE.currentfile.ID);
end


function removesite_callback(a,b,obj)
ind=obj.guihandles.sitelist.Value;
for k=1:length(ind)
    siteID(k)=obj.SE.sites(ind(k)).ID;
end
for k=1:length(ind)
    obj.SE.removeSite(siteID(k));
end
redraw_sitelist(obj);
sv=obj.guihandles.sitelist.Value-1;
sv=min(max(1,sv),length(obj.guihandles.sitelist.String));
% if obj.guihandles.sitelist.Value>length(obj.guihandles.sitelist.String)
obj.guihandles.sitelist.Value=sv;
% end
end

function removecell_callback(a,b,obj);
ind=obj.guihandles.celllist.Value;
cellID=obj.SE.cells(ind).ID;
obj.SE.removeCell(cellID);
redraw_celllist(obj);
redraw_sitelist(obj);
obj.SE.currentcell=obj.SE.cells(1);
% obj.guihandles.sitelist.Value=1;
% obj.guihandles.celllist.Value=1;
end

function angle_callback(data,b,obj)
hax=obj.guihandles.siteax;
sx=(hax.XLim(2)-hax.XLim(1))/4;
alpha=str2double(data.String);
pos=obj.SE.currentsite.pos/1000;
posr=rotatepos([pos(1)-sx pos(2); pos(1)+sx pos(2)],pos,-alpha);
obj.SE.currentsite.annotation.rotationpos=posr;
end



function plotline(obj,posfield,anglehandle)
switch posfield
    case 'line1'
        color='c';
    case 'line2'
        color='m';
    otherwise
        color='b';
end
if nargin<3||isempty(anglehandle)
    anglehandle=[];
end
hax=obj.guihandles.siteax;
alphaimage=obj.SE.currentsite.image.angle;
site=obj.SE.currentsite;
pos=site.annotation.(posfield).pos;
% posfield
% site.annotation.(posfield)
if isa(obj.hlines.(posfield),'imline')
    delete(obj.hlines.(posfield))
end
if sum(pos(:).^2)==0
    obj.hlines.(posfield)=imline(hax);
%         obj.hlines.(posfield).wait
    posin=obj.hlines.(posfield).getPosition;
    roipositon(posin,obj,posfield,anglehandle,hax);
else
%     alphaimage
    posr=rotatepos(pos,site.pos/1000,alphaimage);
%      posr(2)
    obj.hlines.(posfield)=imline(hax,posr);
    roipositon(posr,obj,posfield,anglehandle,hax);

end
obj.hlines.(posfield).setColor(color);

addNewPositionCallback( obj.hlines.(posfield),@(pos) roipositon(pos,obj,posfield,anglehandle,hax));
end

function anglebutton_callback(data,b,obj)
plotline(obj,'rotationpos',obj.guihandles.angle);
end

function roipositon(pos,obj,posfield,outputhandle,outputaxis)
alphaimage=obj.SE.currentsite.image.angle;
angle=pos2angle(pos)+alphaimage;
obj.SE.currentsite.annotation.(posfield).pos=rotatepos(pos,obj.SE.currentsite.pos/1000,-alphaimage);
obj.SE.currentsite.annotation.(posfield).angle=angle;
% pos
len=sqrt(sum((pos(1,:)-pos(2,:)).^2))*1000;
obj.SE.currentsite.annotation.(posfield).length=len;
obj.SE.currentsite.annotation.(posfield).value=len;
% rotatepos(pos,obj.SE.currentsite.pos/1000,-alphaimage)
if nargin>3&&~isempty(outputhandle)&&ishandle(outputhandle)
   
    outputhandle.String=[num2str(angle,'%3.1f') ', ' num2str(len,'%3.0f')];
end

if nargin>4
if ishandle(obj.hlines.sidemarker.(posfield))
    delete(obj.hlines.sidemarker.(posfield))
end
obj.hlines.sidemarker.(posfield)=text(pos(1,1),pos(1,2),'#','Parent',outputaxis,'Color','w','HitTest','off','PickableParts','none','FontSize',15);
end

end


function posr=rotatepos(pos,posroi,alphad)
alpha=alphad/180*pi;
x=pos(:,1);
y=pos(:,2);
x0=x-posroi(1);
y0=y-posroi(2);
posr=zeros(size(pos));
posr(:,1)=x0*cos(alpha)+y0*sin(alpha)+posroi(1);
posr(:,2)=-x0*sin(alpha)+y0*cos(alpha)+posroi(2);

end

function infostruct=getinfostruct(file)
fid=fopen(file);
line=fgetl(fid);
ind=1;
while line>0
    s=eval(['{' line '}']);
    infostruct(ind).text=s{1};
    infostruct(ind).field=s{2};
    ind=ind+1;
    line=fgetl(fid);
end
fclose(fid);
end

function keypress(a,event,obj,flag)
switch event.Key
    case 'rightarrow'
        obj.nextsite(1);
    case 'leftarrow'
        obj.nextsite(-1);
end

end