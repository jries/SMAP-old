classdef SiteExplorer<interfaces.GuiModuleInterface & interfaces.LocDataInterface
    properties
%         sePar
        sites
        numberOfSites=0;
        maxsite=0;
        cells
        numberOfCells=0;
        numberOfFiles=0;
        maxcell=0;
        processors
        files={};
%         currentsiteID=0;
%         currentcellID=0;
%         currentfileID=1;
        currentsite
        currentcell
        currentfile       
        temp

    end
    methods
        function obj=SiteExplorer(varargin)
%             if nargin>0
                obj@interfaces.LocDataInterface(varargin{:});
%             end
            obj.sites=interfaces.SEsites.empty;
            obj.cells=interfaces.SEsites.empty;
            obj.files=interfaces.SEsites.empty;
            obj.temp.siteincell=[];
            obj.temp.cellinfile=[];
        end
        function clear(obj)
            obj.sites=interfaces.SEsites.empty;
            obj.cells=interfaces.SEsites.empty;
            obj.files=interfaces.SEsites.empty;
            obj.numberOfSites=0;
            obj.maxsite=0;
            obj.numberOfCells=0;
            obj.numberOfFiles=0;
            obj.maxcell=0;

        end
%         function attachLocData(obj,locData)
%             attachLocData@interfaces.LocDataInterface(obj,locData);
%             obj.initializeProcessors()
%         end
        function addSites(obj, SEin,filenew,files)
            if isempty(filenew)
                filenew=1:max([SEin.files(:).ID]);
            end
            
            %renumber files
            if ~isempty(SEin.files)
            for k=1:length(SEin.files)
                conversion(SEin.files(k).ID)=filenew(k);
                SEin.files(k).ID=filenew(k);%filenew(SEin.files(k).ID);
                
            end
            %rename files
            for k=1:length(SEin.files)
                SEin.files(k).name=files(k).name;
            end
            
            for k=1:length(SEin.cells)
                SEin.cells(k).info.filenumber=conversion(SEin.cells(k).info.filenumber);
                SEin.cells(k).ID=SEin.cells(k).ID+obj.maxcell;
            end
            for k=1:length(SEin.sites)
                SEin.sites(k).info.filenumber=conversion(SEin.sites(k).info.filenumber);
                SEin.sites(k).info.cell=SEin.sites(k).info.cell+obj.maxcell;
                SEin.sites(k).ID=SEin.sites(k).ID+obj.maxsite;
            end
            
            
            %files
            SEin.numberOfFiles=length(SEin.files);
            for k=1:length(SEin.files)
                obj.files(k+obj.numberOfFiles)=SEin.files(k);
            end
            obj.numberOfFiles=obj.numberOfFiles+SEin.numberOfFiles;
            
%           cells
            SEin.numberOfCells=length(SEin.cells);
            for k=1:length(SEin.cells)
                obj.cells(k+obj.numberOfCells)=SEin.cells(k);
            end
%             obj.cells(obj.numberOfCells+1:obj.numberOfCells+SEin.numberOfCells)=SEin.cells;
            obj.numberOfCells=obj.numberOfCells+SEin.numberOfCells;
%            sites
            SEin.numberOfSites=length(SEin.sites);
            for k=1:length(SEin.sites)
                obj.sites(k+obj.numberOfSites)=SEin.sites(k);
            end
%             obj.sites(obj.numberOfSites+1:obj.numberOfSites+SEin.numberOfSites)=SEin.sites;
            obj.numberOfSites=obj.numberOfSites+SEin.numberOfSites;            
            
            if ~isempty(obj.sites)
            obj.maxsite=max([obj.sites.ID]);
            end
            if ~isempty(obj.cells)
            obj.maxcell=max([obj.cells.ID]);
            end
            
            if ~isempty(obj.cells)
              obj.currentcell=obj.cells(1);
            end
            if ~isempty(obj.files)
              obj.currentfile=obj.files(1);
            end
            if ~isempty(obj.sites)
              obj.currentsite=obj.sites(1);
            end
            
            
%             obj.currentsite=obj.sites(1);
            obj.setIndList;
            end
            
        end
        function setIndList(obj)
            for k=1:length(obj.sites)
                obj.sites(k).indList=k;
            end
        end
        function ID=addSite(obj,site)
            if (isempty(site.ID)||site.ID==0)
                obj.numberOfSites=obj.numberOfSites+1;
                obj.maxsite=obj.maxsite+1;
                site.ID=obj.maxsite;
                site.indList=obj.numberOfSites;
                obj.sites(obj.numberOfSites)=site;
                
            else 
                disp('site already there, ID not 0')
                
            end
            ID=site.ID;
        end
        function removeSite(obj,siteID)
            ind=obj.indexFromID(obj.sites,siteID);
            obj.numberOfSites=obj.numberOfSites-1;
            obj.sites(ind)=[];
        end
        function removeCell(obj,cellID)
            ind=obj.indexFromID(obj.cells,cellID);
            obj.cells(ind)=[];
            si=[obj.sites.info];
            if ~isempty(si)
            indout=[si.cell]==cellID;
            obj.sites(indout)=[];
            obj.numberOfSites=length(obj.sites);
            end
            obj.numberOfCells=length(obj.cells);
            
        end
        function removeFile(obj,fileID)
            indout=obj.indexFromID(obj.files,fileID);
            si=[obj.files.ID];
            indrename=find(si>fileID);
            for ir=indrename
                obj.files(ir).ID=obj.files(ir).ID-1;
            end
            obj.files(indout)=[];
            
            %update sites
            si=[obj.sites.info];
            if ~isempty(si)
                indrename=find([si.filenumber]>fileID);
                indout=find([si.filenumber]==fileID);
                for ir=indrename
                    obj.sites(ir).info.filenumber=obj.sites(ir).info.filenumber-1;
                end
                obj.sites(indout)=[];
                obj.numberOfSites=obj.numberOfSites-sum(indout);

                %update cells
                si=[obj.cells.info];
                indrename=find([si.filenumber]>fileID);
                indout=find([si.filenumber]==fileID);
                for ir=indrename
                    obj.cells(ir).info.filenumber=obj.cells(ir).info.filenumber-1;
                end
                obj.cells(indout)=[];
                obj.numberOfCells=obj.numberOfCells-sum(indout);
            end
        end
        function ID=addCell(obj,cell)
            obj.numberOfCells=obj.numberOfCells+1;
            obj.maxcell=obj.maxcell+1;
            cell.ID=obj.maxcell;
            obj.cells(obj.numberOfCells)=cell;
            ID=cell.ID;
            
%             %boxes
%             cellsize=obj.sePar.Settings.cellfov;
%             if isempty(obj.cellboxes.boxsizenm)||obj.cellboxes.boxsizenm~=cellsize
%                 obj.cellboxes.init(
%             end
            
            
        end
        function addFile(obj,name,ID,info)
            ind=obj.indexFromID(obj.files,ID);
            if isempty(ind)
                obj.numberOfFiles=obj.numberOfFiles+1;
                ind=obj.numberOfFiles;
            end
            
            file=interfaces.SEsites;
            file.name=name;
            file.ID=ID;
            file.info=info;
            obj.files(ind)=file;
        end
        
%         function initializeProcessors(obj)
%             obj.processors.renderer=Renderer(obj.locData);
%             obj.processors.drawer=Drawer(obj.locData);
%             obj.processors.displayer=Displayer(obj.locData);
%         end
        
        function imout=plotsite(obj,site,hax,hbox,haxz)
            if nargin<5
                haxz=[];
            end
            if isnumeric(site)
                ind=obj.indexFromID(obj.sites,site);
                site=obj.sites(ind);
            end
            
            if isempty(site.image)
                
%                  display('draw site')
    %             p1=obj.locData.parameters;
                if length(site.pos)<3
                    site.pos(3)=0;
                end
                p1.sr_pos=site.pos;
                p1.sr_size=ones(2,1)*obj.getPar('se_sitefov')/2;
                p1.sr_pixrec=obj.getPar('se_sitepixelsize');
    %             p1.pixrec=obj.sePar.Settings.sitepixelsize;
                p1.sr_sizeRecPix=round((p1.sr_size*2)/p1.sr_pixrec);
    %             p1.sr_axes=hax;
                p1.sr_axes=-1;
                p1.normalizeFoV=p1.sr_sizeRecPix(1)/obj.getPar('se_sitefov')*obj.getPar('se_siteroi')/2;
                angle=pos2angle(site.annotation.rotationpos.pos);
                if obj.getPar('se_rotate')&&angle~=0
                 p1.rotationangle=angle;
                else 
                 p1.rotationangle=0;
                end

                if obj.getPar('se_imaxcheck_site')
                    p1.imaxtoggle=false;
                    imax=obj.getPar('se_imax_site');

                    for k=1:length(imax)
                        pl{k}.imax_min=imax(k);
                    end
                    for k=length(imax)+1:obj.getPar('numberOfLayers')
                        pl{k}.imax_min=imax(1);
                    end
                else
                    pl=[];
                end

    %              figure(89);
    %             haxz=gca;

                if ~isempty(haxz)
                    p1.sr_size(3)=obj.getPar('se_dz')/2;
                    [site.image, imz]=obj.plotobject(p1,site.info.filenumber,pl);%filenumber
                else
                    site.image=obj.plotobject(p1,site.info.filenumber,pl);%filenumber
                end
                site.image.angle=p1.rotationangle; %remove later? not needed
                if ishandle(haxz)
                    displayimage(imz,haxz);
                    axis(haxz,'equal')
                    set(haxz,'YDir','normal')
                    line(site.pos(1)/1000+[-1 1]*p1.sr_size(1),[site.pos(3) site.pos(3)]/1000,'Color',[1 1 1],'Parent',haxz,'LineWidth',1)
%                     plotbox
                end
            end
            
            if nargin>2&&ishandle(hax)
             displayimage(site.image,hax);
             plotbox(hax,site.pos,obj.getPar('se_siteroi'));
             plotcirc(hax,site.pos,obj.getPar('se_siteroi'));
             line(site.pos(1)/1000+[-.001 .001],site.pos(2)/1000+[0 0],'Color',[1 1 1],'Parent',hax,'LineWidth',3)
             delete(obj.temp.siteincell);
             obj.temp.siteincell=plotbox(hbox,site.pos,obj.getPar('se_sitefov'));
            end
           

            
            imout=site.image;
            
        end
        
        function imout=plotcell(obj,cell,hax,hbox)
             if isempty(cell.image)
                p1.sr_pos=cell.pos;
                p1.sr_size=ones(2,1)*obj.getPar('se_cellfov')/2;
                p1.sr_pixrec=obj.getPar('se_cellpixelsize');
                p1.sr_sizeRecPix=round((p1.sr_size*2)/p1.sr_pixrec);       
                p1.sr_axes=hax;
                p1.sr_axes=-1;
                p1.rotationangle=0;
                p1.normalizeFoV=[];
                
                if obj.getPar('se_imaxcheck_cell')
                    p1.imaxtoggle=false;
                    imax=obj.getPar('se_imax_cell');

                    for k=1:length(imax)
                        pl{k}.imax_min=imax(k);
                    end
                    for k=length(imax)+1:obj.getPar('numberOfLayers')
                        pl{k}.imax_min=imax(1);
                    end
                else
                    pl=[];
                end
                
                
                image=obj.plotobject(p1,cell.info.filenumber,pl);%filenumber
                image.image=single(image.image);
                cell.image=copyfields([],image,{'image','rangex','rangey','parameters','layers','imax'});
             end
            displayimage(cell.image,hax)
            
            
             %plot all site boxes
             if obj.getPar('se_drawboxes')&&~isempty(obj.sites)&&~isempty(obj.sites(1).info)
                 allsites=[obj.sites(:).info];
                 ind=[allsites.cell]==cell.ID;
                 use=getFieldAsVector(obj.sites,'annotation','use');
%                  use(isnan(use))=false;
                if iscell(use)
                    use=ind;
                end
                 
                 plotmanyboxes(hax,obj.sites,ind&use,obj.getPar('se_sitefov'),[1 0 1]);
                 plotmanyboxes(hax,obj.sites,ind & ~use,obj.getPar('se_sitefov'),[1 0 0]);
             end
            
            delete(obj.temp.cellinfile);
            obj.temp.cellinfile=plotbox(hbox,cell.pos,obj.getPar('se_cellfov'));
            
            imout=cell.image;
            
           
            
        end
        
        function plotfile(obj,fileID,hax) %same as cell: store image, temp. check
%             p=obj.getAllParameters;
            ind=obj.indexFromID(obj.files,fileID);
            file=obj.files(ind);
            if isempty(file.image)
                if isfield(file.info,'cam_pixelsize_um')
                pixrec=file.info.cam_pixelsize_um*1000;
                else
                    pixrec=obj.locData.files.file(file.ID).info.cam_pixelsize_um*1000;
                end
                roi=file.info.roi;
                roi([1 3])=roi([1 3])*pixrec(1);
                roi([2 4])=roi([2 4])*pixrec(2);

                % test for SML could be all wrong
                p1.sr_pos=[roi(1)+roi(3)/2 roi(2)+roi(4)/2];
                p1.sr_size=[roi(3) roi(4)]/2; 
                p1.sr_pixrec=mean(pixrec);
                p1.sr_sizeRecPix=round((p1.sr_size*2)/p1.sr_pixrec);
                p1.sr_axes=-1;
                p1.mingausspix=0.8;
                p1.gaussfactor=0.1;
                p1.rotationangle=0;
                p1.normalizeFoV=[];
                image=obj.plotobject(p1,fileID);
                image.image=single(image.image);
                file.image=copyfields([],image,{'image','rangex','rangey','parameters'});
            end
            displayimage(file.image,hax);
            
                         %plot all site boxes
             if obj.getPar('se_drawboxes')&&~isempty(obj.cells)&&~isempty(obj.cells(1).info)
                 allcells=[obj.cells(:).info];
                 ind=[allcells.filenumber]==fileID;
                 plotmanyboxes(hax,obj.cells,ind,obj.getPar('se_cellfov'));
             end
            
        end
        
        function [image,imagez]=plotobject(obj,p,filenumber,pl)
            imagez=[];
           fileind=obj.indexFromID(obj.files,filenumber);
%            p1=obj.locData.parameters;
           numlayers=obj.getPar('numberOfLayers');
           allfields=[renderSMAP drawerSMAP displayerSMAP];
%            allfields=[obj.processors.renderer.inputParameters obj.processors.drawer.inputParameters obj.processors.displayer.inputParameters];
           players=obj.getLayerParameters(1:numlayers,allfields);
           if ~iscell(players)
               players={players};
           end
          plotz=false;
          imax=zeros(numlayers,1);
           for k=1:numlayers
%                 pr=obj.getLayerParameters(k, obj.processors.renderer.inputParameters);   
%                 pd=obj.getLayerParameters(k, obj.processors.drawer.inputParameters); 
%                 pr=copyfields(pr,p);pd=copyfields(pd,p);

                pr=copyfields(players{k},p);
                if nargin>3&&~isempty(pl)&&length(pl)>=k
                    pr=copyfields(pr,pl{k});
                end
                if pr.layercheck
                    pr.ch_filelist.Value=fileind;
                    pr.ch_filelist.selection=pr.ch_filelist.String{fileind};

%                     obj.processors.renderer.setParameters(pr)
%                     obj.processors.drawer.setParameters(pr);
                
                %filter filenumber
%                     groupc=pr.groupcheck;
                    filterold=obj.locData.getFilter(k);
%                     filternew=filterold;
%                     locs=obj.locData.getloc({'filenumber','xnm','ynm'},'grouping',groupc);
                    obj.locData.filter('filenumber',k,'inlist',filenumber)
                    obj.locData.filter('xnm',k,'minmax',[p.sr_pos(1)-p.sr_size(1),p.sr_pos(1)+p.sr_size(1)])
                    obj.locData.filter('ynm',k,'minmax',[p.sr_pos(2)-p.sr_size(2),p.sr_pos(2)+p.sr_size(2)])
                    
                
%                     filternew.filenumber=(locs.filenumber==filenumber);
%                     filternew.xnm=rec.LocalizationFilter.minMaxFilter(locs.xnm,p.sr_pos(1)-p.sr_size(1),p.sr_pos(1)+p.sr_size(1));
%                     filternew.ynm=rec.LocalizationFilter.minMaxFilter(locs.ynm,p.sr_pos(2)-p.sr_size(2),p.sr_pos(2)+p.sr_size(2));
%                 filternew=myrmfield(filternew,'xnm');
%                 filternew=myrmfield(filternew,'ynm');
%                 obj.locData.setFilter(filternew,k)
                
                    
                    
%                     plotz=true;
                   
                    posh=[pr.sr_pos(1) pr.sr_pos(2) pr.sr_size(1)*2 pr.sr_size(2)*2];
                    locz=obj.locData.getloc({'xnm','ynm','znm','locprecnm','locprecznm','PSFxnm','intensity_render','phot','numberInGroup',pr.renderfield.selection},'layer',k,'position',posh);
                    if strcmpi('tiff', pr.rendermode.selection)
                        rawimage=renderSMAP(obj.locData,pr,k);
                    else
                    rawimage=renderSMAP(locz,pr,k);
                    end
%                      rawimage=renderSMAP(obj.locData,pr,k);
                    
                    if ~isempty(locz.znm)&&nargout>1
                        locz.x=locz.xnm;locz.y=locz.znm;
                        locz.sx=locz.locprecnm;locz.sy=locz.locprecznm;
                        prz=pr;
                        prz.normalizeFoV=[];
                        if length(prz.sr_pos)<3
                            prz.sr_pos(2)=0;
                        else
                            
                            prz.sr_pos(2)=prz.sr_pos(3);
                        end
                        prz.sr_size(2)=prz.sr_size(3);
                        rawimagez=renderSMAP(locz,prz,k);
                        layersz(k).images.finalImages=drawerSMAP(rawimagez,prz);
                        plotz=true;
                    end
                    
                    obj.locData.setFilter(filterold,k);
                    layers(k).images.finalImages=drawerSMAP(rawimage,pr);
                    imax(k)=layers(k).images.finalImages.imax;
                    layers(k).images.rawimage=rawimage;
%                      obj.processors.displayer.setParameters(pr);
                    layers(k).images.renderimages=displayerSMAP(layers(k),pr);
                end
           end
%            pd=obj.getLayerParameters(k, obj.processors.displayer.inputParameters); 
%            pd=copyfields(pd,p);

           image=displayerSMAP(layers,pr);
           image.parameters=pr;
           image.parameters.layerparameters=players;
           image.layers=layers;
           image.imax=imax;
           
           if plotz
               imagez=displayerSMAP(layersz,prz);
           end
%            axis equal
            
        end
           
        
        function ind=indexFromID(obj,list,ID)
%             ind=[];
%             if isfield(list,'ID')
            ind=find([list.ID]==ID,1,'first');
%             end
        end
        function cell=getcell(obj,ID)
            ind=obj.indexFromID(obj.cells,ID);
            cell=obj.cells(ind);
        end
        function sites=getsite(obj,ID)
            ind=obj.indexFromID(obj.sites,ID);
            sites=obj.sites(ind);
        end
        function updateSite(obj)
            if obj.currentsite.ID>0
                obj.sites(obj.indexFromID(obj.sites,obj.currentsite.ID))=obj.currentsite;
            end
        end
        function updateCell(obj)
            if obj.currentcell.ID>0
                obj.cells(obj.indexFromID(obj.cells,obj.currentcell.ID))=obj.currentcell;
            end
        end
        
        
        function se=save(obj)
%             se.sePar=obj.sePar;
            se.sites=copy(obj.sites);
            for k=1:length(se.sites)
                se.sites(k).handles=[];
                se.sites(k).image=[];
            end
            se.cells=copy(obj.cells);
            for k=1:length(se.cells)
                se.cells(k).handles=[];
                se.cells(k).image=[];
            end
            
            se.files=copy(obj.files);
            for k=1:length(se.files)
                se.files(k).image=[];
                se.files(k).handles=[];
            end
%             
%             se.files.handles=[];
            se.numberOfSites=obj.numberOfSites;
            se.numberOfCells=obj.numberOfCells;
            se.numberOfFiles=obj.numberOfFiles;
        end
        function load(obj,se)
        end
    end
end

function displayimage(img,hax)
 imagesc(img.rangex,img.rangey,img.image,'Parent',hax,'Pickable','none','HitTest','off')

set(hax,'Xlim',double(img.rangex))
set(hax,'Ylim',double(img.rangey))
set(hax,'YDir','reverse')
hax.HitTest='on';
hax.PickableParts='all';
end

function hg=plotbox(h,pos,roi,color)
if nargin<4
    color=[1 1 1];
end
linewidth=1;

pos=pos/1000;
roi=roi/1000;
x1=pos(1)-roi/2;
x2=pos(1)+roi/2;
y1=pos(2)-roi/2;
y2=pos(2)+roi/2;
hg=hggroup('Parent',h);
line([x1 x1],[y1 y2],'Parent',hg,'Color',color,'LineWidth',linewidth);
line([x1 x2],[y2 y2],'Parent',hg,'Color',color,'LineWidth',linewidth);
line([x2 x2],[y2 y1],'Parent',hg,'Color',color,'LineWidth',linewidth);
line([x2 x1],[y1 y1],'Parent',hg,'Color',color,'LineWidth',linewidth);
end

function hg=plotcirc(h,pos,roi,color)
if nargin<4
    color=[1 1 1]*.5;
end
linewidth=1;

pos=pos/1000;
roi=roi/1000;
x1=pos(1)-roi/2;
x2=pos(1)+roi/2;
y1=pos(2)-roi/2;
y2=pos(2)+roi/2;
hg=rectangle('Parent',h,'Position',[x1,y1,roi,roi],'Curvature',[1,1],'EdgeColor',color,'LineStyle','--');
end

function hg=plotmanyboxes(hax,list,ind,roi,color)
hg=hggroup('Parent',hax);
if nargin<5
color=[1 0 1];
end
for k=1:length(list)
    if ind(k)
        plotbox(hg,list(k).pos,roi,color);
    end
end
end