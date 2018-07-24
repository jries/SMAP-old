classdef AverageSites_kinetochore<interfaces.DialogProcessor&interfaces.SEProcessor
    properties
        imloc
        impos
    end
    methods
        function obj=AverageSites_kinetochore(varargin)        
                obj@interfaces.DialogProcessor(varargin{:});
            obj.showresults=true;
        end
        
        function out=run(obj,p)  
            out=[];
            
%             if p.addfile
%                 obj.locData.addfile(p.name);
%                 initGuiAfterLoad(obj);
%                 obj.SE.processors.preview.updateFilelist;
%             end
%             newfile=obj.locData.files.filenumberEnd;
%             locnew=obj.locData.loc;
            sites=obj.locData.SE.sites;
%             used=false(size(locnew.xnm));
%             x0=nanmedian(locnew.xnm);
%             y0=nanmedian(locnew.ynm);
%             if ~isfield(locnew,'class')
%                 locnew.class=0*locnew.xnm;
%             end
            
            siteline=p.siteline.selection;
            cellline=setdiff(p.siteline.String,siteline);cellline=cellline{1};
            
            window=max(p.size);
            pixelsize=1;
            range=-window:pixelsize:window;
            imc=zeros(length(range)-1,length(range)-1,3);
            imcenters=zeros(length(range)-1,length(range)-1,3);
            
            h=fspecial('gauss',ceil(3*p.filter/pixelsize),p.filter/pixelsize);
            ticc=tic;
            for k=1:length(sites)
                if p.sortselection.Value==1 ||  sites(k).annotation.use
                    [locs1,indloc]=obj.locData.getloc({'xnm','ynm'},'position',sites(k),'layer',1);
                    [locs2,indloc]=obj.locData.getloc({'xnm','ynm'},'position',sites(k),'layer',2);
                    switch p.referencepoint.selection
                        case 'ch1'
                            cp=sites(k).annotation.(siteline).pos(1,:)*1000;
                        case 'ch2'
                            cp=sites(k).annotation.(siteline).pos(2,:)*1000;
                        case 'mean'
                            cp=mean(sites(k).annotation.(siteline).pos*1000,1);
                    end
                    x1=locs1.xnm-cp(1);y1=locs1.ynm-cp(2);
                    [x1r,y1r]=rotcoord(x1,y1,sites(k).annotation.(cellline).angle*pi/180);
                    x2=locs2.xnm-cp(1);y2=locs2.ynm-cp(2);
                    [x2r,y2r]=rotcoord(x2,y2,sites(k).annotation.(cellline).angle*pi/180);
                    
                    
                     xl1=sites(k).annotation.(siteline).pos(1,1)*1000-cp(1);yl1=sites(k).annotation.(siteline).pos(1,2)*1000-cp(2);
                    [xlr1,ylr1]=rotcoord(xl1,yl1,sites(k).annotation.(cellline).angle*pi/180);
                     xl2=sites(k).annotation.(siteline).pos(2,1)*1000-cp(1);yl2=sites(k).annotation.(siteline).pos(2,2)*1000-cp(2);
                    [xlr2,ylr2]=rotcoord(xl2,yl2,sites(k).annotation.(cellline).angle*pi/180);
%                     
                    im1=histcounts2(x1r,y1r,range,range);
                    im2=histcounts2(x2r,y2r,range,range);
                      imc(:,:,1)=imc(:,:,1)+im1;imc(:,:,2)=imc(:,:,2)+im2;
 
                    im1centers=histcounts2(xlr1,ylr1,range,range);
                    im2centers=histcounts2(xlr2,ylr2,range,range);
                    imcenters(:,:,1)=imcenters(:,:,1)+im1centers;imcenters(:,:,2)=imcenters(:,:,2)+im2centers;
%                       imc(:,:,3)=imc(:,:,3)+iml;
                      
                      
%                     locnew.xnm(indloc)=locs.xnm-sites(k).pos(1);
%                     locnew.ynm(indloc)=locs.ynm-sites(k).pos(2);
%     %                 figure(88)
%     %                 plot(locnew.xnm(indloc),locnew.ynm(indloc),'+')
%                     locnew.filenumber(indloc)=newfile;
%                     locnew.class(indloc)=sites(k).ID;
%                     used=used|indloc;
                end
                if toc(ticc)>1
                    ticc=tic;
                    obj.status(['average site: ' num2str(k) ' of ' num2str(length(sites))]); drawnow
                end
            end
            prof=squeeze(sum(imc,2));
            
            if p.filterc
            imc=imfilter(imc,h);
             imcenters=imfilter(imcenters,h);
            end
            ax=obj.initaxis('localizations');
            imagesc(ax,range,range,imc/max(imc(:)));
            ax=obj.initaxis('profile');
            plot(ax,range(1:end-1),prof)
            ax=obj.initaxis('centers');
            
            switch p.referencepoint.selection
               case 'ch1'
                    imcenterss=imcenters(:,:,2);
               case 'ch2'
                    imcenterss=imcenters(:,:,1);
               case 'mean'
                   imcenterss=imcenters;
            end
            
           
            imagesc(ax,range,range,imcenterss/max(imcenterss(:)))
            obj.imloc=imc/max(imc(:));
            obj.impos=imcenterss/max(imcenterss(:));
            
%             afsd
%             used=used& ~( isnan(locnew.xnm) | isnan(locnew.ynm) | isnan(locnew.class));
%             
%             fn=fieldnames(locnew);
%               for k=1:length(fn)
%                    locnew.(fn{k})=locnew.(fn{k})(used);
%               end
%              
%              locc=obj.locData.copy;
%              locc.loc=locnew;
%              locc.regroup;
%              locc.filter;
% 
%             if p.addfile
%                 locnew.xnm=locnew.xnm+x0;
%                 locnew.ynm=locnew.ynm+y0;
%             
%                 for k=1:length(fn)
%                     obj.locData.addloc(fn{k},locnew.(fn{k}))
%                 end
%                 obj.locData.regroup;
%                 obj.locData.filter;
%             end
%            
%             
%             
%             %try: add empty file, there put averaged sites
%             %for every site: loc.xnm-site.pos(1)+xpossite
%             
%             loc1=locc.getloc({'xnm','ynm'},'layer',1,'position','all','filenumber',newfile);
%             loc2=locc.getloc({'xnm','ynm'},'layer',2,'position','all','filenumber',newfile);
%             %do some statistics:
%            
%             
%             
%             ax=obj.initaxis('scatter');
%             plot(loc1.xnm,loc1.ynm,'.')
%             hold on
%             plot(loc2.xnm,loc2.ynm,'.')
%             hold off
%             
%              [phi1,r1]=cart2pol(loc1.xnm,loc1.ynm);
%              [phi2,r2]=cart2pol(loc2.xnm,loc2.ynm);
%             ax=obj.initaxis('radial distribution');
%             dr=5;
% %             rr=dr/2:dr:max(max(r1),max(r2));
%             %proper concentration: dA=pi*(ro^2-ri^2)
%             rr=0:dr:max(max(r1),max(r2))+dr;
%             
%             dA=pi*(rr(2:end)+rr(1:end-1)).*(rr(2:end)-rr(1:end-1));
%             h1=histcounts(r1,[rr]) ;
%              h2=histcounts(r2,[rr]) ;
%               rrp=rr(1:end-1)+dr/2;
%             plot(rrp,h1,rrp,h2)
%             xlabel('r')
%             ylabel('counts')
%             
%             ax=obj.initaxis('radial concentration');
% %             plot(rr,h1./rr/2/pi,rr,h2./rr/2/pi)
%            
%              plot(rrp,h1./dA,rrp,h2./dA)
%             xlabel('r')
%             ylabel('concentration')
        end
        function pard=guidef(obj)
            pard=guidef(obj);
        end
        function savebutton_callback(obj,a,b)
            fl=obj.getPar('mainfile');
            [path, file]=fileparts(fl);
            [file, path]=uiputfile([path file '.tif']);
            
            if file
                fn1=[path  strrep(file,'.tif','_locs.tif')];
                fn2=[path  strrep(file,'.tif','_centers.tif')];
                imwrite(obj.imloc,fn1);
                imwrite(obj.impos,fn2);
            end
            
        end
    end
end




function pard=guidef(obj)

pard.t1.object=struct('String','average kinetochores','Style','text');
pard.t1.position=[1,1];
pard.t1.Width=1;

pard.sitelinet.object=struct('String','Site line','Style','text');
pard.sitelinet.position=[2,1];
pard.sitelinet.Width=1;

pard.siteline.object=struct('String',{{'line1','line2'}},'Style','popupmenu');
pard.siteline.position=[2,2];
pard.siteline.Width=1;

pard.referencepointt.object=struct('String','ref point','Style','text');
pard.referencepointt.position=[2,3];
pard.referencepointt.Width=1;

pard.referencepoint.object=struct('String',{{'ch1','ch2','mean'}},'Style','popupmenu');
pard.referencepoint.position=[2,4];
pard.referencepoint.Width=1;

pard.sizet.object=struct('String','size x y (nm)','Style','text');
pard.sizet.position=[3,1];
pard.sizet.Width=1;

pard.size.object=struct('String','100 50','Style','edit');
pard.size.position=[3,2];
pard.size.Width=1;

pard.post.object=struct('String','all or annotated:use','Style','text');
pard.post.position=[4,1];
pard.post.Width=1;


pard.sortselection.object=struct('String',{{'all','use'}},'Style','popupmenu');
pard.sortselection.position=[4,2];
pard.sortselection.Width=1;

pard.filterc.object=struct('String','Filter (size nm):','Style','checkbox');
pard.filterc.position=[5,1];
pard.filterc.Width=1;

pard.filter.object=struct('String','10','Style','edit');
pard.filter.position=[5,2];
pard.filter.Width=1;

pard.savebutton.object=struct('String','save images','Style','pushbutton','Callback',@obj.savebutton_callback);
pard.savebutton.position=[5,4];
pard.savebutton.Width=1;

% pard.namet.object=struct('String','name','Style','text');
% pard.namet.position=[5,1];
% pard.namet.Width=1;
% 
% pard.name.object=struct('String','average','Style','edit');
% pard.name.position=[5,2];
% pard.name.Width=1;
% 
% pard.addfile.object=struct('String','add average as new data set','Style','checkbox');
% pard.addfile.position=[5,3];
% pard.addfile.Width=2;

pard.plugininfo.type='ROI_Analyze';


end