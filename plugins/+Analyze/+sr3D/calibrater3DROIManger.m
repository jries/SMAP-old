classdef calibrater3DROIManger<interfaces.DialogProcessor
    properties
        SXY
    end
    methods
        function obj=calibrater3DROIManger(varargin)   
            obj@interfaces.DialogProcessor(varargin{:}) ;
             obj.showresults=true;
        end
        function initGui(obj)
%             beaddistribution_callback(0,0,obj)           
        end
        function out=run(obj,p)
            out=[];
            
            [path,file]=fileparts(obj.getPar('lastSMLFile'));
            file=strrep(file,'_sml','_3Dcal');
%             [file,path]=uiputfile([path filesep file]);
            roisize=p.roisize;
            halfroisize=(roisize-1)/2;
            roifac=1.5;
            roisizebig=(roifac*roisize);
            halfroisizebig=round((roisizebig-1)/2); %room for shifting
            roisizebig=2*halfroisizebig+1;
%             framerange=p.framerange;
            framerange=1:max(obj.locData.loc.frame);
            
            
            se=obj.locData.SE;
            sites=se.sites;
            sitefilenumbers=getFieldAsVector(se.sites,'info','filenumber');
            usesites=getFieldAsVector(se.sites,'annotation','use');
            
            %f0 for all beads: 
            ax=obj.initaxis('find f0');
            for k=length(sites):-1:1
                locs=obj.locData.getloc({'frame','PSFxnm','PSFynm','phot'},'position',sites(k));
                f0(k)=getf0site(locs);
            end
            
            %if on glass
            fzero=round(median(f0(~isnan(f0))));
            
            allrois=zeros(roisizebig,roisizebig,length(framerange),sum(usesites));
            files=se.files;
            induse=1;
            k=1; %for loop
                filename=files(k).info.basefile;
                roi=files(k).info.roi;
                campix=files(k).info.pixsize;
                il=getimageloader(obj,filename);
                imstackadu=il.getmanyimages([],'mat');
                
                imstack=(double(imstackadu)-il.metadata.offset)*il.metadata.pix2phot;
                sim=size(imstack);
                thisfile=find(sitefilenumbers==k);
                for s=1:length(thisfile)
                    sitenumber=thisfile(s);
                    if usesites(sitenumber)
                        pos=round(sites(sitenumber).pos(1:2)/campix/1000-roi(1:2));
                        if pos(1)>halfroisizebig&&pos(1)<sim(1)-halfroisizebig&&pos(2)>halfroisizebig&&pos(2)<sim(2)-halfroisizebig
                            allrois(:,:,:,induse)=imstack(pos(2)-halfroisizebig:pos(2)+halfroisizebig,pos(1)-halfroisizebig:pos(1)+halfroisizebig,framerange);
                            induse=induse+1;
                        end
                    end
                end
                
            ax=obj.initaxis('CC PSF');
            [corrPSF,shiftedstack]=registerPSFsCorr(allrois,fzero);
            drawnow;
            corrPSFhd=zeros(size(corrPSF,1)*2,size(corrPSF,2)*2,size(corrPSF,3));
            for k=1:size(corrPSF,3)
                corrPSFhd(:,:,k)=imresize(corrPSF(:,:,k),2);
            end
            s_size=p.roisize;
            tic
            [np_psf,coeff]=generate_psf_to_spline_YLJ(corrPSFhd,s_size,fzero);
            toc
%             bspline
            np_psf = np_psf/max(np_psf(:));
            for i = 1:size(np_psf,3)
            np_psf(:,:,i) = np_psf(:,:,i)/sum(sum(np_psf(:,:,i)));
            end

%             tic
%             b3_0=bsarray(double(np_psf),'lambda',p.smoothingfactor);
%             toc
            tic
            b3_0=bsarray(double(corrPSFhd),'lambda',p.smoothingfactor);
            toc
            
%             out.coeff=coeff;
            bspline=b3_0;
            save([path filesep file],'bspline','coeff')
%             b3_0_1=bsarray(double(np_psf),'lambda',0.1);%with smoothing factor 0.1
            

            %upsampling
            %spline generation
            ax=obj.initaxis('PSFz')
            ftest=fzero-15;
            zpall=squeeze(shiftedstack(halfroisizebig+1,halfroisizebig+1,:,:));
            xpall=squeeze(shiftedstack(:,halfroisizebig+1,ftest,:));
            
            for k=1:size(zpall,2)
                zpall(:,k)=zpall(:,k)/max(zpall(:,k));
                xpall(:,k)=xpall(:,k)/max(xpall(:,k));
            end
            
            zprofile=squeeze(corrPSF(halfroisizebig+1,halfroisizebig+1,:));
            zprofile=zprofile/max(zprofile);
            
            
            xprofile=squeeze(corrPSF(:,halfroisizebig+1,ftest));
            xprofile=xprofile/max(xprofile);
            
            [XX,YY,ZZ]=meshgrid(1:b3_0.dataSize(1),1:b3_0.dataSize(2),framerange);
             psfbs = interp3_0(b3_0,XX,YY,ZZ);
             
            zbs=squeeze(psfbs(2*halfroisizebig+1,2*halfroisizebig+1,:));
            zbs=zbs/max(zbs);
            xbs=squeeze(psfbs(:,2*halfroisizebig+1,ftest));
            xbs=xbs/max(xbs);
            
           
%             [XX,YY,ZZ]=meshgrid(1:roisizebig*2,halfroisizebig*2+1,ftest);
%             psfbs = interp3_0(b3_0,XX,YY,ZZ);
%             xbs=squeeze(psfbs(:));
%             xbs=xbs/max(xbs);
%             fbs=fzero-13:1:fzero+12;
            
            hold off
            plot(framerange,zpall,'c')
            hold on
            plot(framerange,zprofile,'k+')
            plot(framerange,zbs,'ro')
            
            
            xrange=-halfroisizebig:halfroisizebig;
            
             ax=obj.initaxis('PSFx')
            hold off
            plot(xrange,xpall,'c')
            hold on
            plot(xrange,xprofile,'k+-')
            xrangebig=-halfroisizebig-.25:0.5:halfroisizebig+.5;
            plot(xrangebig,xbs,'ro')
            
        end
        function pard=guidef(obj)
            pard=guidef(obj);
        end
    end
end

function f0=getf0site(locs)
[zas,zn]=stackas2z(locs.PSFxnm,locs.PSFynm,locs.frame,locs.phot,true);
drawnow
f0=zas;
end

function il=getimageloader(obj,filename)
try
    il=imageloaderAll(filename,[],obj.P); %still exist?
catch
    maindir=obj.getGlobalSetting('DataDirectory');
    filename=strrep(filename,'\','/');
    d=dir(maindir);
    alldir={d([d.isdir]).name};
    ind=[1 strfind(filename,'/')];
    for k=1:length(ind)-1
        thisf=filename(ind(k)+1:ind(k+1)-1);
        if any(contains(alldir,thisf))&&~isempty(thisf)


            filename=[maindir '/' filename(ind(k)+1:end)];
            break
        end
    end
    try
        il=imageloaderAll(filename,[],obj.P); %still exist?
    catch
        [~,~,ext]=fileparts(filename);
        if isempty(ext)
            filename=[filename '.tif'];
        end
        [f,path]=uigetfile(filename);
        if f
            il=imageloaderAll([path f],[],obj.P); 
        else
            il=[];
        end
    end
    %look for file in main directory
end
end


function pard=guidef(obj)
pard.dzt.object=struct('Style','text','String','dz (nm)'); 
pard.dzt.position=[1,1];
pard.dzt.Width=.5;
pard.dz.object=struct('Style','edit','String','50'); 
pard.dz.position=[1,1.5];
pard.dz.Width=.5;

pard.beaddistribution.object=struct('String',{{'Beads on Glass','Beads in Gel'}},'Style','popupmenu','Callback',{{@beaddistribution_callback,obj}});
pard.beaddistribution.position=[2,1];
pard.beaddistribution.Width=2;


pard.roisizet.object=struct('Style','text','String','Roi (pix)'); 
pard.roisizet.position=[3,1];
pard.roisizet.Width=.7;
pard.roisize.object=struct('Style','edit','String','13'); 
pard.roisize.position=[3,1.7];
pard.roisize.Width=.3;

pard.framewindowt.object=struct('Style','text','String','frame window size'); 
pard.framewindowt.position=[3,2];
pard.framewindowt.Width=1;
pard.framewindow.object=struct('Style','edit','String','15'); 
pard.framewindow.position=[3,3];
pard.framewindow.Width=.5;

pard.smoothingfactort.object=struct('Style','text','String','Smoothing factor Bspline'); 
pard.smoothingfactort.position=[3,3.5];
pard.smoothingfactort.Width=1;
pard.smoothingfactor.object=struct('Style','edit','String','.1'); 
pard.smoothingfactor.position=[3,4.5];
pard.smoothingfactor.Width=.5;


% tp=3.1;tmin=3.6;td=3.9;tmax=4.2;
% w=0.3;
% wp=0.5;
% pard.tt1.object=struct('String','grid','Style','text');
% pard.tt1.position=[2,tp];
% pard.tt1.Width=wp;
% pard.mint.object=struct('String','min','Style','text');
% pard.mint.position=[2,tmin];
% pard.mint.Width=w;
% pard.dxt.object=struct('String','delta','Style','text');
% pard.dxt.position=[2,td];
% pard.dxt.Width=w;
% pard.maxt.object=struct('String','max','Style','text');
% pard.maxt.position=[2,tmax];
% pard.maxt.Width=w;
% 
% pard.Xt.object=struct('String','X (pix)','Style','text');
% pard.Xt.position=[3,tp];
% pard.Xt.Width=wp;
% 
% pard.Xmin.object=struct('String','0','Style','edit');
% pard.Xmin.position=[3,tmin];
% pard.Xmin.Width=w;
% pard.Xd.object=struct('String','128','Style','edit');
% pard.Xd.position=[3,td];
% pard.Xd.Width=w;
% pard.Xmax.object=struct('String','512','Style','edit');
% pard.Xmax.position=[3,tmax];
% pard.Xmax.Width=w;
% 
% pard.Yt.object=struct('String','Y (pix)','Style','text');
% pard.Yt.position=[4,tp];
% pard.Yt.Width=wp;
% 
% pard.Ymin.object=struct('String','0','Style','edit');
% pard.Ymin.position=[4,tmin];
% pard.Ymin.Width=w;
% pard.Yd.object=struct('String','128','Style','edit');
% pard.Yd.position=[4,td];
% pard.Yd.Width=w;
% pard.Ymax.object=struct('String','512','Style','edit');
% pard.Ymax.position=[4,tmax];
% pard.Ymax.Width=w;
% 
% pard.Zt.object=struct('String','Z vals (nm)','Style','text');
% pard.Zt.position=[5,tp];
% pard.Zt.Width=wp;
% pard.Zval.object=struct('String',' 0:1000:3000','Style','edit');
% pard.Zval.position=[5,tmin];
% pard.Zval.Width=5-tp-wp;
% 
% % pard.Zmin.object=struct('String','0','Style','edit');
% % pard.Zmin.position=[5,tmin];
% % pard.Zmin.Width=w;
% % pard.Zd.object=struct('String','1000','Style','edit');
% % pard.Zd.position=[5,td];
% % pard.Zd.Width=w;
% % pard.Zmax.object=struct('String','3000','Style','edit');
% % pard.Zmax.position=[5,tmax];
% % pard.Zmax.Width=w;
% 
% pard.splinet.object=struct('String','Spline','Style','text');
% pard.splinet.position=[6,1];
% pard.splinet.Width=1;
% pard.splinepart.object=struct('String','smoothing par','Style','text');
% pard.splinepart.position=[6,2];
% pard.splinepart.Width=1;
% pard.splinepar.object=struct('String','0.95','Style','edit');
% pard.splinepar.position=[6,3];
% pard.splinepar.Width=1;
% 
% pard.fitzt.object=struct('String','fit z MLE:','Style','text');
% pard.fitzt.position=[7,1];
% pard.fitzt.Width=1;
% pard.fitzranget.object=struct('String','z range (nm)','Style','text');
% pard.fitzranget.position=[7,2];
% pard.fitzranget.Width=1;
% pard.fitzrange.object=struct('String','-500 500','Style','edit');
% pard.fitzrange.position=[7,3];
% pard.fitzrange.Width=1;
% pard.fitzB0.object=struct('String','B0=0','Style','checkbox','Value',0);
% pard.fitzB0.position=[7,4];
% pard.fitzB0.Width=1;
% 
% pard.calculateZSxSy.object=struct('String','Z(Sx,Sy,X,Y,Z)','Style','checkbox','Value',0);
% pard.calculateZSxSy.position=[8,1];
% pard.calculateZSxSy.Width=1.;
% 
% pard.St.object=struct('String','S (pix)','Style','text');
% pard.St.position=[8,tp];
% pard.St.Width=wp;
% pard.Smin.object=struct('String','0','Style','edit');
% pard.Smin.position=[8,tmin];
% pard.Smin.Width=w;
% pard.Sd.object=struct('String','.03','Style','edit');
% pard.Sd.position=[8,td];
% pard.Sd.Width=w;
% pard.Smax.object=struct('String','4','Style','edit');
% pard.Smax.position=[8,tmax];
% pard.Smax.Width=w;
% 
% pard.framerangecombinet.object=struct('String','zrange beyond interval (nm)','Style','text');
% pard.framerangecombinet.position=[3,1];
% pard.framerangecombinet.Width=1;
% pard.framerangecombine.object=struct('String','100','Style','edit');
% pard.framerangecombine.position=[3,2];
% pard.framerangecombine.Width=.5;
% 
% pard.ztoframet.object=struct('String','set z to frame','Style','checkbox');
% pard.ztoframet.position=[3,1];
% pard.ztoframet.Width=1;
% pard.ztoframe.object=struct('String','21','Style','edit');
% pard.ztoframe.position=[3,2];
% pard.ztoframe.Width=.5;
% 
% 
% pard.zrangeuset.object=struct('String','zrange for fit  (nm)','Style','text');
% pard.zrangeuset.position=[4,1];
% pard.zrangeuset.Width=1;
% pard.zrangeuse.object=struct('String','-800 800','Style','edit');
% pard.zrangeuse.position=[4,2];
% pard.zrangeuse.Width=1;
% 
% pard.getcoords.object=struct('String','select','Style','pushbutton','Callback',{{@getcoords,obj}});
% pard.getcoords.position=[4,4.5];
% pard.getcoords.Width=.5;
% pard.getcoords.Height=2;
% 
% pard.save.object=struct('String','save','Style','pushbutton','Callback',{{@save_SXY,obj}});
% pard.save.position=[9,1];
% pard.save.Width=1;
% pard.save.Height=1;
% 
% pard.correctbeadsizec.object=struct('String','Correc for bead size (nm)','Style','checkbox','Value',0);
% pard.correctbeadsizec.position=[1,3];
% pard.correctbeadsizec.Width=1.5;
% 
% pard.beadsize.object=struct('String','100','Style','edit');
% pard.beadsize.position=[1,4.5];
% pard.beadsize.Width=.5;
% 
% % pard.savebutton.object=struct('String','save','Style','pushbutton');
% % pard.savebutton.position=[7,3];
pard.inputParameters={'cam_pixelsize_nm'};
pard.plugininfo.type='ProcessorPlugin';


end
