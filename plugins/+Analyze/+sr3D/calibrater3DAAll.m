classdef calibrater3DAAll<interfaces.DialogProcessor
    properties
        SXY
    end
    methods
        function obj=calibrater3DAAll(varargin)   
            obj@interfaces.DialogProcessor(varargin{:}) ;
             obj.showresults=true;
             obj.guiPar.FieldHeight=obj.guiPar.FieldHeight-1;obj.guiPar.Vrim=obj.guiPar.Vrim-20;
        end
        function out=run(obj,p)
            out=[];
            %get beads. Either parse sites or segment new
 %                 beads.loc.xnm,ynm,frame,psfxnm,psfynm
%                 beads.filenumber .stack  .f0 .pos 
            if p.beadsource.Value==1 %ROImanager
                beads=sites2beads(obj,p);
            else % segment new
                beads=segmentb(obj,p);
            end
            p.EMon=obj.locData.files.file(1).info.EMon;
            p.fminmax=[min(obj.locData.loc.frame) max(obj.locData.loc.frame)];
           axall=getaxes(p);

            % get f0 for beads
            for k=1:length(beads)
                beads(k).loc.PSFxpix=beads(k).loc.PSFxnm/p.cam_pixelsize_nm;
                 beads(k).loc.PSFypix=beads(k).loc.PSFynm/p.cam_pixelsize_nm;
                [beads(k).f0,beads(k).psfx0,beads(k).psfy0]=getf0site(beads(k).loc,p);
            end
            badind=isnan([beads(:).f0]);
            beads(badind)=[];
            
            p.fglass=myquantile([beads(:).f0],0.03);
            %get image stacks if needed
            if p.fitbsplinec || p.fitcsplinec 
                beads=getimagestacks(obj,p,beads);      
            end
            
            %if iteratively: for k=1:2
            if p.refine 
                redo=2;
            else
                redo=1;
            end
            
            p=getranges(p);
            indgood=true(1,length(beads));
            
            for iter=1:redo  
                beads=beads(indgood);
                indgood=true(1,length(beads));
            %if spatial analysis: group according to X,Y,Z 
            %get clean curves
            xbead=zeros(length(beads),1);
            ybead=zeros(length(beads),1);
            for k=length(beads):-1:1
                xbead(k)=beads(k).pos(1);
                ybead(k)=beads(k).pos(2);
            end
            for X=1:length(p.Xrange)-1
                for Y=1:length(p.Yrange)-1
                    indh=xbead>p.Xrange(X)-p.xyoverlap&xbead<p.Xrange(X+1)+p.xyoverlap&ybead>p.Yrange(Y)-p.xyoverlap&ybead<p.Yrange(Y+1)+p.xyoverlap;
                    if sum(indh)==0
                        continue
                    end
                    beadh=beads(indh);
                    % get calibration and used beads
                    
                    [curvecal,indgoodc]=getcurvecal(beadh,p,X,Y,axall);
                    if p.fitbsplinec|| p.fitcsplinec
                        [stackcal,indgoods]=getstackcal(beadh,p,X,Y,axall);
                        for Z=1:length(curvecal)
                            allcal(Z)=copyfields(stackcal(Z),curvecal(Z));
                        end
                    else
                        allcal=curvecal;
                        indgoods=indgoodc;
                    end
    
                    
                    indhf=find(indh);
                    indgood(indhf)=indgood(indhf)&indgoods&indgoodc;
%                     indgood=indgood&indgoods&indgoodc;
                    SXY(X,Y,:)=allcal;
                end
            end       
            end
            plotcurves(obj,SXY,axall,p)
            %save
            [path,file]=fileparts(obj.getPar('lastSMLFile'));
            file=strrep(file,'_sml','_3Dcal');
            save([path filesep file],'SXY')
        end
        
        function initGui(obj)
            setvisible(0,0,obj)
%             beaddistribution_callback(0,0,obj)           
        end
 
        function pard=guidef(obj)
            pard=guidef(obj);
        end
    end
end

%%
function bead=segmentb(obj,p)

minframes=(p.zrangeuse(2)-p.zrangeuse(1))/ p.dz/2;

locDatacopy=obj.locData.copy;
locDatacopy.regroup(150,minframes);

beadind=find(locDatacopy.grouploc.numberInGroup>minframes);
xb=locDatacopy.grouploc.xnm(beadind);
yb=locDatacopy.grouploc.ynm(beadind);
filenumber=locDatacopy.grouploc.filenumber(beadind);
winsize=250;
for k=length(beadind):-1:1
%     beadhere=(mywithin(locDatacopy.loc.xnm,[xb(k)-winsize/2 winsize],locDatacopy.loc.ynm,[yb(k)-winsize/2 winsize]));
%     beadfile=locDatacopy.loc.filenumber==filenumber(k);
%     beadhere=beadhere&beadfile;
    bead(k).loc=locDatacopy.getloc({'xnm','ynm','PSFxnm','PSFynm','frame','phot'},'filenumber',filenumber(k),'Position',[xb(k) yb(k) winsize]);
%     bead(k).loc=copystructReduce(locDatacopy.loc,beadhere,{'xnm','ynm','PSFxnm,','PSFynm','frame'});
    bead(k).filenumber=filenumber(k);
    bead(k).pos=[xb(k) yb(k)];
    
end
end

function bead=sites2beads(obj,p)

se=obj.locData.SE;
sites=se.sites;
for k=length(sites):-1:1
    bead(k).loc=obj.locData.getloc({'xnm','ynm','frame','PSFxnm','PSFynm','phot'},'grouping','ungrouped','position',sites(k));
    bead(k).filenumber=sites(k).info.filenumber;
    bead(k).pos=sites(k).pos;   
end
use=logical(getFieldAsVector(sites,'annotation','use'));
bead=bead(use);
end

function [f0,PSFx0,PSFy0]=getf0site(locs,p)
p.ploton=0;
[zas,zn]=stackas2z(locs.PSFxpix,locs.PSFypix,locs.frame,locs.phot,p);
if isnan(zas)
    PSFx0=NaN;
    PSFy0=NaN;
else
    ind=find(locs.frame<=zas,1,'last');
PSFx0=locs.PSFxpix(ind);
PSFy0=locs.PSFypix(ind);
end
% drawnow
f0=zas;
end

function beads=getimagestacks(obj,p,beads)

roisize=p.roisize;
% halfroisize=(roisize-1)/2;
roifac=1.5;
roisizebig=(roifac*roisize);
halfroisizebig=round((roisizebig-1)/2); %room for shifting
roisizebig=2*halfroisizebig+1;
storeframes=p.roiframes+500/p.dz;
halfstoreframes=round((storeframes-1)/2);
storeframes=halfstoreframes*2+1;
f0=[beads(:).f0];
f0r=round(f0);

fminmax=[min(obj.locData.loc.frame) max(obj.locData.loc.frame)];
if p.beaddistribution.Value==2
    fzero=halfstoreframes+1;
else%if on glass
    fzero=round(median(f0(~isnan(f0))));
end

framerange0=max(fminmax(1),fzero-halfstoreframes):min(fzero+halfstoreframes,fminmax(2));

files=obj.locData.files.file;
indbad=false(length(beads),1);
 for k=1:length(files) 
    if isfield(files(k).info,'imagefile')
        filename=files(k).info.imagefile;
    else
        filename=files(k).info.basefile;
    end
    roi=files(k).info.roi;
    campix=files(k).info.pixsize;
    il=getimageloader(obj,filename);
    imstackadu=il.getmanyimages([],'mat');
                
    imstack=(double(imstackadu)-il.metadata.offset)*il.metadata.pix2phot;
    sim=size(imstack);
    thisfile=find([beads(:).filenumber]==k);
%                 f=figure(88);
%             hold off    
%             imagesc(imstack(:,:,40));
            
    for s=length(thisfile):-1:1
        beadnumber=thisfile(s);
        pos=round(beads(beadnumber).pos(1:2)/campix/1000-roi(1:2));
        if pos(1)>halfroisizebig&&pos(1)<sim(1)-halfroisizebig&&pos(2)>halfroisizebig&&pos(2)<sim(2)-halfroisizebig
            if p.beaddistribution.Value==1 %on glass
                frange=framerange0;
            else
                frange=max(fminmax(1),f0r(beadnumber)-halfstoreframes):min(fminmax(2),f0r(beadnumber)+halfstoreframes);
            end
            smallf=imstack(pos(2)-halfroisizebig:pos(2)+halfroisizebig,pos(1)-halfroisizebig:pos(1)+halfroisizebig,frange);
            stack=zeros(roisizebig,roisizebig,storeframes);
            % XXXXX gel: or always: center fzero or similar.
            if frange(1)==fminmax(1)
                stack(:,:,end-size(smallf,3)+1:end)=smallf;
            else
                stack(:,:,1:size(smallf,3))=smallf;
            end
            beads(beadnumber).stack.image=stack;
            beads(beadnumber).stack.framerange=frange;
            
%             figure(88)
% hold on
% plot(pos(1),pos(2),'wo')
%             figure(89)
%             imagesc(smallf(:,:,12));
%             title(pos)
            
            
%             imageslicer(stack,'Parent',f)
            
        else
            indbad(beadnumber)=true;
        end
                    
    end
 end
beads(indbad)=[];
end

function p=getranges(p)

if p.spatialcalibration
    
    Xrange=(p.Xmin:p.Xd:p.Xmax).*p.cam_pixelsize_nm;
    Yrange=(p.Ymin:p.Yd:p.Ymax).*p.cam_pixelsize_nm;

    Zrange=p.Zval;
    if ischar(Zrange)
        Zrange=str2num(Zrange)
    end
    if length(Xrange)==1
        Xrange(2)=p.Xmax*p.cam_pixelsize_nm;
    end
    if length(Yrange)==1
        Yrange(2)=p.Ymax*p.cam_pixelsize_nm;
    end
    if length(Zrange)==1
        Zrange(2)=p.Zmax;
    end
    p.Zrange=Zrange;
    p.Xrange=Xrange;
    p.Yrange=Yrange;
    if p.beaddistribution.Value==1 %Glass
        p.Zrange=[-inf inf];
    end
    
else
    p.Zrange=[-inf inf];
    p.Xrange=[-inf inf];
    p.Yrange=[-inf inf];   
end
end

function axall=getaxes(p)
            ax=initaxis(p.resultstabgroup,'beads');
            axall.htbeads=uitabgroup(ax.Parent);
            axall.axbeads=maketgax(axall.htbeads,'scatter');   
            
            if p.fitbsplinec||p.fitcsplinec
                ax=initaxis(p.resultstabgroup,'spline fit');
                axall.allsplines=uitabgroup(ax.Parent);
                ax=maketgax(axall.allsplines,'scatter');
%                 ax=initaxis(axall.allsplines,'scatter');
                axall.hspline_scatter=uitabgroup(ax.Parent);
                
                ax=maketgax(axall.allsplines,'PSF');
                axall.hspline_psf=uitabgroup(ax.Parent);
                
                ax=maketgax(axall.allsplines,'Overlay');
                axall.hspline_overlay=uitabgroup(ax.Parent);
                
                 ax=maketgax(axall.allsplines,'PSFz');
                axall.hspline_psfz=uitabgroup(ax.Parent);
                
                ax=maketgax(axall.allsplines,'PSFx');
                axall.hspline_psfx=uitabgroup(ax.Parent);
                
                ax=maketgax(axall.allsplines,'validate');
                axall.hspline_validate=uitabgroup(ax.Parent);
            end  
                
            %init validation and summary axes
            ax=initaxis(p.resultstabgroup,'splines');
            axall.htsplines=uitabgroup(ax.Parent);
            axall.axsplines=maketgax(axall.htsplines,'summary'); 

            if p.fitzsxsyc
                ax=initaxis(p.resultstabgroup,'sx^2-sy^2');
                axall.htsx2sy2=uitabgroup(ax.Parent);
                axall.axsxsys=maketgax(axall.htsx2sy2,'summary');  
                axall.axsxs22=maketgax(axall.htsx2sy2,'validation'); 
            end
            if p.fitzc
                ax=initaxis(p.resultstabgroup,'z fit');
                axall.htzfit=uitabgroup(ax.Parent);
                axall.axzfits=maketgax(axall.htzfit,'summary');
            end
            if p.calculateZSxSy
                ax=initaxis(p.resultstabgroup,'Z(Sx,Sy)');
                axall.hzsx=uitabgroup(ax.Parent);
 
                axall.axzlut=maketgax(axall.hzsx,'validation'); 
            end
              
end

function [SXY,indgoodc]=getcurvecal(beadsh,p,X,Y,axall)
%z-dependent? 
zc=p.spatialcalibration && p.zcalc &p.beaddistribution.Value==2;

% hold off
indgoodc=true(1,length(beadsh));
for Z=1:length(p.Zrange)-1
    for B=length(beadsh):-1:1
        beadz0=(beadsh(B).f0-p.fglass)*p.dz;
        
        if p.alignz||p.beaddistribution.Value==2
             beadz=(beadsh(B).loc.frame*p.dz)-beadz0;
        else 
            beadz=(beadsh(B).loc.frame)*p.dz;
        end
        zglass=p.fglass*p.dz;
        if zc
            if p.zfilter.Value==1 % f0

                if beadz0>p.Zrange(Z)-p.framerangecombine && beadz0<p.Zrange(Z+1)+p.framerangecombine
                   sx=beadsh(B).loc.PSFxpix; 
                   sy=beadsh(B).loc.PSFypix; 
                   z=beadz-zglass;     
                   phot=beadsh(B).loc.phot;  
                else
                    sx=[];sy=[];z=[];
                end             
            else
                zs=(beadsh(B).loc.frame-p.fglass)*p.dz;
                goodz=zs>p.Zrange(Z)-p.framerangecombine & zs<p.Zrange(Z+1)+p.framerangecombine;
               sx=beadsh(B).loc.PSFxpix(goodz); 
               sy=beadsh(B).loc.PSFypix(goodz); 
               z=beadz(goodz)-zglass;   
               phot=beadsh(B).loc.phot(goodz);
                
            end
        else
           sx=beadsh(B).loc.PSFxpix; 
           sy=beadsh(B).loc.PSFypix; 
           z=beadz-zglass; phot=beadsh(B).loc.phot;  
        end
        inzr=z>=p.zrangeuse(1)&z<=p.zrangeuse(2);
        curves(B).sx=double(sx(inzr));
        curves(B).sy=double(sy(inzr));
        curves(B).z=double(z(inzr));
        curves(B).phot=double(phot(inzr));
        curves(B).xpos=beadsh(B).pos(1);
        curves(B).ypos=beadsh(B).pos(2);
%         plot(curves(B).z,curves(B).sx,curves(B).z,curves(B).sy)
%         hold on
    end
    
    %get calibrations
    bh=curves;
    tgt=[num2str(X) num2str(Y) num2str(Z)];
    ax=maketgax(axall.htsplines,tgt);                     
     [spline,indg]=getcleanspline(curves,p);
     indgoodc=indgoodc&indg;
    bh=bh(indg);
    
    SXY(Z)=getsxyinit(p,X,Y,Z);
    SXY(Z).spline=spline;
    if p.fitzsxsyc
        ax=maketgax(axall.htsx2sy2,tgt);
        SXY(Z).Sx2_Sy2=cal_Sx2_Sy2(bh,p);
    end
    if p.fitzc
        ax=maketgax(axall.htzfit,tgt);
        SXY(Z).fitzpar=cal_fitzpar(bh,p);
    end
    if p.calculateZSxSy
        ax=maketgax(axall.hzsx,tgt);
        SXY(Z).splineLUT=cal_splineLUT(SXY(Z).spline,p);
    end
    SXY(Z).curve=curves;

end
end



       function out=run(obj,p)
%             out=[];
%             
%             
%             locsall=obj.locData.getloc({'frame','xnm','ynm','PSFxnm','PSFynm','filenumber','phot'},'position','all','layer',1,'removeFilter','filenumber','grouping','ungrouped');
% %             locsall=obj.locData.getloc('frame','xnm','ynm','PSFxnm','PSFynm','filenumber','phot');
%             if isempty(locsall.PSFynm)
%                 error('no PSFy found')
%             end
%             if p.correctbeadsizec
%                 Rbead=p.beadsize/2;
%                 x=-Rbead:Rbead;
%                 y=sqrt(Rbead.^2-x.^2);
%                 sigmaB=sqrt(sum(y.*x.^2)/sum(y));
%                 locsall.PSFxnm=sqrt(locsall.PSFxnm.^2-sigmaB.^2);
%                 locsall.PSFynm=sqrt(locsall.PSFynm.^2-sigmaB.^2);
%             end
%             
%             locsall.PSFxpix=locsall.PSFxnm/p.cam_pixelsize_nm;
%             locsall.PSFypix=locsall.PSFynm/p.cam_pixelsize_nm;
%             locsall.zfnm=locsall.frame*p.dz;
%             
%             ax=initaxis(p.resultstabgroup,'found beads');
%             htg=uitabgroup(ax.Parent);
%                          axbeadss=maketgax(htg,'all');
%                          axbeadss.NextPlot='add';
%             numfiles=max(locsall.filenumber);
%             maxd=p.cam_pixelsize_nm*2;
%             locsall.beadnum=zeros(size(locsall.xnm));
%             for k=1:numfiles
%                 ax=maketgax(htg,num2str(k));
%                 indf=locsall.filenumber==k;
%                 beadlocs(k)=getBeadLocs(locsall.xnm(indf),locsall.ynm(indf),p);
%                 [beadnum,numlocs]=associatelocs(beadlocs(k).x,beadlocs(k).y,locsall.xnm(indf),locsall.ynm(indf),maxd);
%                 beadn=beadnum>0;
%                 indff=find(indf);
%                 locsall.beadnum(indff(beadn))=beadnum(beadn)+max(locsall.beadnum);
%                 plot(axbeadss,beadlocs(k).x,beadlocs(k).y,'o')
%             end
%             
%             %make bead structure
%             for k=max(locsall.beadnum):-1:1
%                 thisbead=(locsall.beadnum==k);
%                 bead(k).PSFxpix=double(locsall.PSFxpix(thisbead));
%                 bead(k).PSFypix=double(locsall.PSFypix(thisbead));
%                 bead(k).phot=double(locsall.phot(thisbead));
%                 bead(k).xpos=median(double(locsall.xnm(thisbead)));
%                 bead(k).ypos=median(double(locsall.ynm(thisbead)));
%                 bead(k).filenumber=double(locsall.filenumber(find(thisbead,1)));
%                 bead(k).numlocs=sum(thisbead);
%                 bead(k).frame=double(locsall.frame(thisbead));
%                 bead(k).zfnm=double(locsall.zfnm(thisbead));
%                 bead(k).beadnum=double(locsall.beadnum(thisbead));
%                 
%                 %derived properties
%                 bead(k).ztrue=stackas2z(bead(k).PSFxpix,bead(k).PSFypix,bead(k).zfnm,bead(k).phot,0);
%                 indz0=find(bead(k).zfnm>bead(k).ztrue,1);
%                 rangez0=indz0-15:indz0+15;rangez0(rangez0<1)=1;rangez0(rangez0>length(bead(k).zfnm))=length(bead(k).zfnm);
%                 if isempty(rangez0)
%                     bead(k).minSx=NaN;
%                     bead(k).minSy=NaN;
%                 else
%                     bead(k).minSx=min(bead(k).PSFxpix(rangez0));
%                     bead(k).minSy=min(bead(k).PSFypix(rangez0));                    
%                 end
% %                 ztrue(k);
%                 bead(k).I0=max(bead(k).phot);
%             end
%             
%             % if on glass: correct position based on average, not taking
%             % into account position. This might create errors. change
%             % later?
%             zrangeall=[min(locsall.frame) max(locsall.frame)].*p.dz;
%             if p.beaddistribution.Value==1 %glass
%                 f1ind=find([bead(:).filenumber]==1);
%                 zf1=robustMean([bead(f1ind).ztrue]);
%                 if p.ztoframet
%                     zpos=p.ztoframe.*p.dz;
%                 else
%                     zpos=0;
%                 end
% %                 zshift=dzh+zpos;
%                 zrangeall=zrangeall-zpos;
%                 for k=1:max([bead(:).filenumber])
%                         bi=find([bead(:).filenumber]==k);
%                         dzh=robustMean([bead(bi).ztrue])-zf1;
%                         for b=1:length(bi)
% %                             dzh+zpos
%                             bead(bi(b)).ztrue=bead(bi(b)).ztrue-dzh-zpos;
%                             bead(bi(b)).zfnm=bead(bi(b)).zfnm-dzh-zpos;
%                             bead(bi(b)).zrangeall=zrangeall;
%                         end
%                 end
%                 p.ztruepos=zf1-zpos; 
%                 p.Zval=[0 0];
%             else
%                 %make absolute with respect to coverslip. Either from GUI
%                 %or from finding smallest ztrue. Assume: approx same stack
%                 %z position (rest taken care of by ztrue)
% %                 zshift=0;
%                 ztrue=[bead(:).ztrue];
%                 zGlass=quantile(ztrue,0.02);
%                 for k=1:length(bead)
%                     bead(k).ztrue=bead(k).ztrue-zGlass;
%                     bead(k).zfnm=bead(k).zfnm-zGlass;
%                     bead(k).zrangeall=zrangeall-zGlass;
%                 end
%                 p.ztruepos=0;
%             end
%             
%             
%             
%             %sort beads according to X,Y
%             
%             Xrange=(p.Xmin:p.Xd:p.Xmax).*p.cam_pixelsize_nm;
%             Yrange=(p.Ymin:p.Yd:p.Ymax).*p.cam_pixelsize_nm;
%             
%             Zrange=p.Zval;
%             if ischar(Zrange)
%                 Zrange=str2num(Zrange)
%             end
%             if length(Xrange)==1
%                 Xrange(2)=p.Xmax*p.cam_pixelsize_nm;
%             end
%             if length(Yrange)==1
%                 Yrange(2)=p.Ymax*p.cam_pixelsize_nm;
%             end
%             if length(Zrange)==1
%                 Zrange(2)=p.Zmax;
%             end
%             p.Zrange=Zrange;
%             
%             %init validation and summary axes
%             ax=initaxis(p.resultstabgroup,'splines');
%             htsplines=uitabgroup(ax.Parent);
%             ax=initaxis(p.resultstabgroup,'sx^2-sy^2');
%             htsx2sy2=uitabgroup(ax.Parent);
%             ax=initaxis(p.resultstabgroup,'z fit');
%             htzfit=uitabgroup(ax.Parent);
%             ax=initaxis(p.resultstabgroup,'Z(Sx,Sy)');
%             hzsx=uitabgroup(ax.Parent);
%             axzfits=maketgax(htzfit,'summary');  
%             axsxsys=maketgax(htsx2sy2,'summary');  
%             axsxs22=maketgax(htsx2sy2,'validation');  
%             axsplines=maketgax(htsplines,'summary'); 
%             axzlut=maketgax(hzsx,'validation'); 
%             
%             %get clean curves
%             for X=1:length(Xrange)-1
%                 for Y=1:length(Yrange)-1
%                     indh=[bead(:).xpos]>Xrange(X)&[bead(:).xpos]<Xrange(X+1)&[bead(:).ypos]>Yrange(Y)&[bead(:).ypos]<Yrange(Y+1);
%                     if sum(indh)==0
%                         continue
%                     end
%                     bh=cleanupbeads(bead(indh),p);
%                     curves=getcurves(bh,p);
%                     for Z=1:length(curves)
%                         CXY(X,Y,Z)=curves(Z);
%                     end
%                 end
%             end           
%             scurves=size(CXY);
%             if length(scurves)==2
%                 scurves(3)=1;
%             end
%             
% %             get fits
%             for X=1:scurves(1)%length(Xrange)-1
%                 for Y=1:scurves(2)%length(Yrange)-1
%                     for Z=1:scurves(3)
%                         bh=CXY{X,Y,Z};
%                         tgt=[num2str(X) num2str(Y) num2str(Z)];
%                         ax=maketgax(htsplines,tgt);                     
%                         [SXY(X,Y,Z).spline,indg]=getcleanspline(bh,p);
%                         bh=bh(indg);
%                         ax=maketgax(htsx2sy2,tgt);
%                         SXY(X,Y,Z).Sx2_Sy2=cal_Sx2_Sy2(bh,p);
%                         ax=maketgax(htzfit,tgt);
%                         SXY(X,Y,Z).fitzpar=cal_fitzpar(bh,p);
%                         if p.calculateZSxSy
%                             ax=maketgax(hzsx,tgt);
%                             SXY(X,Y,Z).splineLUT=cal_splineLUT(SXY(X,Y,Z).spline,p);
%                         else
%                             SXY(X,Y,Z).splineLUT=[];
%                         end
%                         SXY(X,Y,Z).curve=CXY{X,Y,Z};
%                         SXY(X,Y,Z).Xrangeall=Xrange;
%                         SXY(X,Y,Z).Yrangeall=Yrange;
%                         SXY(X,Y,Z).Zrangeall=Zrange;
%                         SXY(X,Y,Z).posind=[X,Y,Z];                        
%                         SXY(X,Y,Z).Xrange=[Xrange(X), Xrange(X+1)];
%                         SXY(X,Y,Z).Yrange=[Yrange(Y) ,Yrange(Y+1)];
%                         SXY(X,Y,Z).Zrange=[Zrange(Z), Zrange(Z+1)];
%                         SXY(X,Y,Z).Zoffset=p.ztruepos+mean(SXY(X,Y,Z).Zrange);
%                         SXY(X,Y,Z).legend=[num2str(X) num2str(Y) num2str(Z)];
%                     end
%                 end
%             end
%             
       end
       function plotcurves(obj,SXY,axall,p)
            %plot results
            axes(axall.axbeads)
            sp=SXY(:);
            legends={SXY(:).legend};
            for k=1:numel(sp)
                if ~isempty(sp(k).curve)
                    xpos=vertcat(sp(k).curve(:).xpos);
                    ypos=vertcat(sp(k).curve(:).ypos);
                plot(xpos,ypos,'k.');hold on;
                text(mean(sp(k).Xrange),mean(sp(k).Yrange),sp(k).legend)
                end
            end
            
            for k=1:length(p.Xrange)
                line([p.Xrange(k),p.Xrange(k)],[p.Yrange(1) p.Yrange(end)])
            end
            for k=1:length(p.Yrange)
                line([p.Xrange(1),p.Xrange(end)],[p.Yrange(k) p.Yrange(k)])
            end
            
            
            sp=[SXY(:).spline];
%             zt=zrangeall(1):0.01:zrangeall(2);
%             z0=zf1-zpos;
%             z0=zshift;
            z0=0;
            zt=z0+p.zrangeuse(1):0.01:z0+p.zrangeuse(2);

            linecolors=lines(length(sp));
            axes(axall.axsplines)
            hold off;
            
            for k=1:numel(sp)
                pl2(k)=plot(axall.axsplines,zt,sp(k).x(zt),'Color',linecolors(k,:));
                hold on;
                plot(axall.axsplines,zt,sp(k).y(zt),'Color',linecolors(k,:));
                
            end
            ylim([0 5])
            legend(pl2,legends);
            
            if ~isempty(SXY(1).Sx2_Sy2)
                sp={SXY(:).Sx2_Sy2};
                s=-30:0.5:30;
                axes(axall.axsxsys)
                hold off;
                for k=1:numel(sp)
                    if ~isempty(sp{k})
                    plot(axall.axsxsys,sp{k}(s),s);hold on;
                    end
                end
                xlabel(axall.axsxsys,'z')
                ylabel(axall.axsxsys,'Sx^2-Sy^2')
                xlim(axall.axsxsys,p.zrangeuse)
                legend(legends);
            end
            
            if ~isempty(SXY(1).fitzpar)
                sp={SXY(:).fitzpar};


                axes(axall.axzfits)
                hold off
                xlabel(axall.axzfits,'z')
                ylabel(axall.axzfits,'Sx, Sy')
                ylim(axall.axzfits,[0 5])

                hold off;
    %             nleg={};
                pl=[];
                for k=1:numel(sp)
                    if ~isempty(sp{k})
                        [sx,sy]=getsxfromzfitpar(zt,sp{k},z0); 
                        pl(k)=plot(axall.axzfits,zt,sx,'Color',linecolors(k,:));
                         hold on
                        plot(axall.axzfits,zt,sy,'Color',linecolors(k,:))

    %                     nleg{k}=num2str(k);
                    end
                end

                legend(pl,legends);
            end
            %cross check and validate with bead positions
%             axes(axzlut);
            
%             hold off
            zt=p.zrangeuse(1):p.dz:p.zrangeuse(2);
            
            for k=1:numel(SXY)
                sxa=vertcat(SXY(k).curve(:).sx);
                sya=vertcat(SXY(k).curve(:).sy);
                za=vertcat(SXY(k).curve(:).z);
                
                if ~isempty(SXY(k).splineLUT)
                    zha=zfromSXSYLut(SXY(k).splineLUT,sxa,sya);
                    dzba=bindata(za,za-zha,zt,'mean');
                    dzs=bindata(za,za-zha,zt,'std');
                    plot(axall.axzlut,za,za-zha,'.','MarkerSize',2)
                    axall.axzlut.NextPlot='add';
                    plot(axall.axzlut,zt,dzba,'k',zt,dzba+dzs,'k',zt,dzba-dzs,'k')
                end
                    zha2=zfromSx2_Ss2(SXY(k).Sx2_Sy2,sxa,sya);
                     plot(axall.axsxs22,za,za-zha2,'.','MarkerSize',2)
                     axall.axsxs22.NextPlot='add';
                    dzba2=bindata(za,za-zha2,zt,'mean');
                    dzs2=bindata(za,za-zha2,zt,'std');
                    plot(axall.axsxs22,zt,dzba2,'k',zt,dzba2+dzs2,'k',zt,dzba2-dzs2,'k')
            end
            plot(axall.axzlut,p.zrangeuse,[0 0],'k');
            plot(axall.axsxs22,p.zrangeuse,[0 0],'k');
            yrange=[-100 100];
            ylim(axall.axzlut,yrange)
            xlim(axall.axzlut,p.zrangeuse)
            ylim(axall.axsxs22,yrange)
            xlim(axall.axsxs22,p.zrangeuse)     
            obj.SXY=SXY;  
       end

function fitzpar=cal_fitzpar(b,p)
 z=vertcat(b(:).z);
Sx=vertcat(b(:).sx);
Sy=vertcat(b(:).sy);
% ztrue=robustMean([b(:).ztrue]);
% zrange=spline.maxmaxrange;
% zrange=zrange+[300 -300];

% midpoint=robustMean([b(:).ztrue]);
% ztruepos=0;
zrange=[p.fitzrange(1) p.fitzrange(2)];
ax=gca;hold off;
fitzpar=getzfitpar(Sx,Sy,z,zrange,0,~p.fitzB0,ax);

end

function splineLUT=cal_splineLUT(spline,p)
srange=p.Smin:p.Sd:p.Smax;
splineLUT=gets2z(spline,srange);
% splineLUT=[];
end

function sxp=cal_Sx2_Sy2(b,p)


z=vertcat(b(:).z);
Sx=vertcat(b(:).sx);
Sy=vertcat(b(:).sy);
% figure(88)
hold off
sxp=fitsx2sy2(Sx,Sy,z,p.zrangeuse,gca);
% sxp.ztruepos=p.ztruepos;
xlim([p.zrangeuse])
end

function [s,indgood2]=getcleanspline(curves,p)

% zm=robustMean([b(:).ztrue]);


 za=vertcat(curves(:).z);
Sxa=vertcat(curves(:).sx);
Sya=vertcat(curves(:).sy);
% indz=za>zm+p.zrangeuse(1)&za<zm+p.zrangeuse(2);
indz=za>p.zrangeuse(1)&za<p.zrangeuse(2);
z=za(indz);
Sx=Sxa(indz);
Sy=Sya(indz);

% Sxs=smoothn({z,Sx});
splinex=fit(z,Sx,'poly6','Robust','LAR','Normalize','on');
spliney=fit(z,Sy,'poly6','Robust','LAR','Normalize','on');

hold off
for k=length(curves):-1:1
    w=(curves(k).phot);
%     w=1;
    zh=curves(k).z;
    indzh=zh>p.zrangeuse(1)&zh<p.zrangeuse(2);
    w=w.*indzh;
    errh=(curves(k).sx-splinex(zh)).^2.*w+(curves(k).sy-spliney(zh)).^2.*w;
    err(k)=sqrt(sum(errh)/sum(w));
    errh2=(curves(k).sx-splinex(zh)).*w./curves(k).sx;
    errh3=(curves(k).sy-spliney(zh)).*w./curves(k).sy;
    err2(k)=abs(sum(errh2)/sum(w));
    err3(k)=abs(sum(errh3)/sum(w));
%     err(k)=sqrt(robustMean(errh));
%     plot(b(k).zfnm,b(k).PSFxpix,'y.')
%     hold on
%     plot(b(k).zfnm,b(k).PSFypix,'y.')
%     plot(b(k).zfnm,w/max(w))
    
end

%  
[em,es]=robustMean(err);
% [em2,es2]=robustMean(err2);
% [em3,es3]=robustMean(err3);
if isnan(es), es=em; end
% if isnan(es2), es2=em2; end
% if isnan(es3), es3=em3; end
indgood2=err<em+2*es & err2+err3<.15;
if sum(indgood2)==0
    indgood2=true(length(curves),1);
end
zg=vertcat(curves(indgood2).z);
indz=zg>p.zrangeuse(1)&zg<p.zrangeuse(2);
zg=zg(indz);
sxg=vertcat(curves(indgood2).sx);
syg=vertcat(curves(indgood2).sy);
sxg=sxg(indz);
syg=syg(indz);

splinex2=getspline(sxg,zg,1./(abs(sxg-splinex(zg))+.1),p.splinepar);
spliney2=getspline(syg,zg,1./(abs(syg-spliney(zg))+.1),p.splinepar);
% splinex2=getspline(sxg,zg,1,p.splinepar);
% spliney2=getspline(syg,zg,1,p.splinepar);
% zt=b(1).zrangeall(1):0.01:b(1).zrangeall(2);
zt=min(zg):0.01:max(zg);

for k=1:length(curves)
    if indgood2(k)
        t='r.';
    else
        t='g.';
    end
    plot(curves(k).z,curves(k).sx,t)
    hold on
    plot(curves(k).z,curves(k).sy,t)
    
end

plot(zt,splinex(zt),'b')
plot(zt,spliney(zt),'b')
plot(zt,splinex2(zt),'k')
plot(zt,spliney2(zt),'k')

xlim([zt(1) zt(end)])
drawnow
s.x=splinex2;s.y=spliney2;
s.zrange=[zt(1) zt(end)];


zr=s.zrange(1):1:s.zrange(2);
midp=round(length(zr)/8);


[~,ind1x]=max(s.x(zr(1:midp)));
[~,ind1y]=max(s.y(zr(1:midp)));
[~,ind2x]=max(s.x(zr(midp:end)));
[~,ind2y]=max(s.y(zr(midp:end)));

z1=max(zr(ind1x),zr(ind1y));
z2=min(zr(ind2x+midp-1),zr(ind2y+midp-1));

s.maxmaxrange=[z1 z2];
end

function spline=getspline(S,z,w,p)
if nargin<4
    p=0.96;
end
% 
[zs,zind]=sort(z);Ss=S(zind);ws=w(zind);
spline=fit(zs,Ss,'smoothingspline','Weights',ws,'Normalize','on','SmoothingParam',p); 
% spline=fit(zs,Ss,'smoothingspline','Weights',ws,'Normalize','on'); 
% spline=fit(zs,Ss,'smoothingspline','Weights',ws,'SmoothingParam',0.99); %weigh smaller more
% [spline,a,b]=fit(zs,Ss,'smoothingspline','Weights',ws); %weigh smaller more
% ds=Ss-spline(zs);
% spline2=fit(zs,Ss,'smoothingspline','Weights',1./ds.^2); %weigh smaller more
end

function calibrateAstig3D(locs,p)
global zt sxf syf
sxt=double(locs.PSFxnm)/p.cam_pixelsize_nm;syt=double(locs.PSFynm)/p.cam_pixelsize_nm;framet=double(locs.frame);
zt=framet*p.dz;

B0=double(~p.B0);


recgui.initaxis(p.resultstabgroup,'select range');
plot(framet,sxt,'ro')
hold on
plot(framet,syt,'bo')
hold off
title('select range (two points)')
[indi,y]=ginput(2);

% indi=round(indum/p.dz*1000);

%rangef=round(max(1,min(indi))):round(min(max(indi),length(sxt)));
rangef=round(max(min(indi))):round(min(max(indi)));

range=find(framet>=rangef(1),1,'first'):find(framet<=rangef(end),1,'last');

sx=sxt(range);sy=syt(range);frame=framet(range);
z=frame*p.dz;



[~,ix]=min(sx);
[~,iy]=min(sy);
% 
midpreal=(z(ix)+z(iy))/2;


indframez0=find(p.framez0<=framet,1,'first');
midp=zt(indframez0);

z=z-midp;zt=zt-midp;
% recgui.initaxis(p.resultstabgroup,'fitted points');

plot(zt,sxt,'r.')
hold on
plot(zt,syt,'b.')
plot(z,sx,'ro')
plot(z,sy,'bo')

mpf=midpreal-midp;
% parx= [d sx0 sy0 Ax Ay Bx By g mp]
startp=[    0.3    1.0    1.0000  0   0        0         0  0.307   -mpf];

% fitp=lsqcurvefit(@sbothfromsigma,startp,[z z],[sx sy])
fitp=lsqnonlin(@sbothfromsigmaerr,startp,[],[],[],[z z],[sx sy],B0);

% fitfunction=@(par,x)sbothfromsigmafit(par,x);
% fit2=fit([z; z],[sx; sy],fitfunction);

sxf=sigmafromz(fitp([1 2 4 6 8 9]),zt,B0);
plot(zt,sxf,'k')
fpy=fitp([1 3 5 7 8 9]);
fpy(5)=-fpy(5);
syf=sigmafromz(fpy,zt,B0);
plot(zt,syf,'k')
hold off

axis tight

% 
% figure(2)
% plot(sx,sy,'.')
% hold on
% plot(sxf,syf,'r')
% plot(syf,sxf,'c')
% hold off


%zpar=[sigma0x,Ax,Ay,Bx,By,gamma,d,sigma0y)
outforfit=real(fitp([2 4 5 6 7 8 1 3]));
button = questdlg('Is the fit good?','3D astigmatism bead calibration','Refit','Save','Cancel','Refit') ;
switch button
    case 'Refit'
        calibrateAstig3D(locs,p)
    case 'Save'
        fn=[p.fileout(1:end-4) '_3DAcal.mat']
        [f,p]=uiputfile(fn);
        if f
            save([p f],'outforfit')
        end
end
end

function s=sigmafromz(par,z,B0)
global gamma
par=real(par);
% parx= [d sx0 Ax Bx g mp]
s0=par(2);d=par(1);A=par(3);B=par(4)*B0;g=par(5);mp=par(6);

% s=s0*sqrt(1+(z-g+mp).^2/d^2);
s=s0*sqrt(1+(z-g+mp).^2/d^2+A*(z-g+mp).^3/d^3+B*(z-g+mp).^4/d^4);
s=real(s);
end

function s=sbothfromsigma(par,z,B0)
% parx= [d sx0 sy0 Ax Ay Bx By g mp]
par=real(par);
px=par([1 2 4 6 8 9]);
py=par([1 3 5 7 8 9]);
 py(5)=-py(5);
zh=z(:,1);
s=real([sigmafromz(px,zh,B0) sigmafromz(py,zh,B0)]);
end

function err=sbothfromsigmaerr(par,z,sx,B0)
sf=sbothfromsigma(par,z,B0);
err=sf-sx;

% stderr=std(err);
err=err./sqrt(abs(err));
end

function getcoords(a,b,obj)
locsall=obj.locData.getloc({'frame','xnm','ynm','PSFxnm','PSFynm','filenumber','phot'},'position','all','layer',1,'removeFilter','filenumber');
p=obj.getAllParameters;
x=locsall.xnm/p.cam_pixelsize_nm;
y=locsall.ynm/p.cam_pixelsize_nm;
img=myhist2(x,y,1,1,[0 512.5],[0 512.5]);
f1=figure;
imagesc(img');
h=imrect;
pos=wait(h);% [x y wx wy]
delete(f1);
answ=inputdlg({'number of rows','number of columns'},'set tiles', 1,{'2','1'});
if isempty(answ)
    return;
end
nx=str2double(answ{1});
ny=str2double(answ{2});

% pos=pos([2 1 4 3]);
p.Xmin=round(pos(1)); p.Ymin=round(pos(2)); p.Xmax=round(pos(1)+pos(3)); p.Ymax=round(pos(2)+pos(4));
p.Xd=floor(pos(3)/nx);
p.Yd=floor(pos(4)/ny);
obj.setGuiParameters(p)

end

function setvisible(a,b,obj)
p=obj.getAllParameters;
% gel={'Zt','Zval'};
glass={'ztoframe','ztoframet','alignz'};
gelpatial={'framerangecombinet','framerangecombine','Zt','Zval','zfilter','zfiltert'};
% p=obj.getSingleGuiParameter('beaddistribution');
switch p.beaddistribution.Value 
    case 1%glass
%         setvis(obj,gel,'off')
        setvis(obj,glass,'on')
        setvis(obj,gelpatial,'off')
        setvis(obj,'zcalc','off')
    case 2 %gel
%         setvis(obj,gel,'on')
        setvis(obj,glass,'off')
        if p.spatialcalibration
            setvis(obj,'zcalc','on')
            if p.zcalc
                setvis(obj,gelpatial,'on')
            else
                setvis(obj,gelpatial,'off')
            end
        else
            setvis(obj,'zcalc','off')
            setvis(obj,gelpatial,'off')
        end
end



% pard.framerangecombinet.object=struct('String','z beyond interval (nm)','Style','text');


scfields={'Ymax','Yd','Ymin','Yt','Xmax','Xd','Xmin','Xt','maxt','dxt','mint','tt1','getcoords','xyoverlap','xyoverlapt'};
if p.spatialcalibration
    setvis(obj,scfields,'on')
else
    setvis(obj,scfields,'off')
end

splinef={'framewindow','alignzt','smoothingfactor','smoothingfactort','roiframes','roiframest','roisize','roisizet'};
if p.fitbsplinec||p.fitcsplinec
    setvis(obj,splinef,'on')
else
    setvis(obj,splinef,'off')
end

gaussz={'fitzB0','fitzrange','fitzranget'};
if p.fitzc
    setvis(obj,gaussz,'on')
else
    setvis(obj,gaussz,'off')
end


spl={'splinepart','splinepar'};
if p.fitzsxsyc||p.calculateZSxSy
    setvis(obj,spl,'on')
else
    setvis(obj,spl,'off')
end

spl={'Smax','Sd','Smin','St'};
if p.calculateZSxSy
    setvis(obj,spl,'on')
else
    setvis(obj,spl,'off')
end

end

function setvis(obj,fields,state)
if ~iscell(fields)
    fields={fields};
end
for k=1:length(fields)
    obj.guihandles.(fields{k}).Visible=state;
end
end

function pard=guidef(obj)
tp=3.6;tmin=4.1;td=4.4;tmax=4.7;
w=0.3;
wp=0.5;
wcb=1.;

pard.beadsource.object=struct('String',{{'RoiManager','Segment'}},'Style','popupmenu');
pard.beadsource.position=[1,1];
pard.beadsource.Width=1;

pard.dzt.object=struct('Style','text','String','dz (nm)'); 
pard.dzt.position=[2,1];
pard.dzt.Width=.4;
pard.dz.object=struct('Style','edit','String','50'); 
pard.dz.position=[2,1.4];
pard.dz.Width=.35;

pard.refine.object=struct('Style','checkbox','String','refine iteratively','Value',0); 
pard.refine.position=[2,2];
pard.refine.Width=1;

pard.beaddistribution.object=struct('String',{{'Glass','Gel'}},'Style','popupmenu','Callback',{{@setvisible,obj}});
pard.beaddistribution.position=[1,2];
pard.beaddistribution.Width=.75;



pard.ztoframet.object=struct('String','z0 (frame)','Style','checkbox');
pard.ztoframet.position=[6,tp];
pard.ztoframet.Width=.75;
pard.ztoframe.object=struct('String','21','Style','edit');
pard.ztoframe.position=[6,tp+.75];
pard.ztoframe.Width=.25;

pard.alignz.object=struct('Style','checkbox','String','Align in z'); 
pard.alignz.position=[7,tp];
pard.alignz.Width=1.;


pard.zrangeuset.object=struct('String','zrange (nm)','Style','text');
pard.zrangeuset.position=[3,1.5];
pard.zrangeuset.Width=.65;
pard.zrangeuse.object=struct('String','-800 800','Style','edit');
pard.zrangeuse.position=[3,2.15];
pard.zrangeuse.Width=.65;


pard.spatialcalibration.object=struct('Style','checkbox','String','Spatial calibration','Value',0,'Callback',{{@setvisible,obj}}); 
pard.spatialcalibration.position=[1,tp];
pard.spatialcalibration.Width=wcb+.2;

pard.getcoords.object=struct('String','select','Style','pushbutton','Callback',{{@getcoords,obj}});
pard.getcoords.position=[1,tp+wcb];
pard.getcoords.Width=3*w+wp-wcb;
% pard.getcoords.Height=1;

pard.tt1.object=struct('String','grid','Style','text');
pard.tt1.position=[2,tp];
pard.tt1.Width=wp;
pard.mint.object=struct('String','min','Style','text');
pard.mint.position=[2,tmin];
pard.mint.Width=w;
pard.dxt.object=struct('String','delta','Style','text');
pard.dxt.position=[2,td];
pard.dxt.Width=w;
pard.maxt.object=struct('String','max','Style','text');
pard.maxt.position=[2,tmax];
pard.maxt.Width=w;

pard.Xt.object=struct('String','X (pix)','Style','text');
pard.Xt.position=[3,tp];
pard.Xt.Width=wp;

pard.Xmin.object=struct('String','0','Style','edit');
pard.Xmin.position=[3,tmin];
pard.Xmin.Width=w;
pard.Xd.object=struct('String','512','Style','edit');
pard.Xd.position=[3,td];
pard.Xd.Width=w;
pard.Xmax.object=struct('String','512','Style','edit');
pard.Xmax.position=[3,tmax];
pard.Xmax.Width=w;

pard.Yt.object=struct('String','Y (pix)','Style','text');
pard.Yt.position=[4,tp];
pard.Yt.Width=wp;

pard.Ymin.object=struct('String','0','Style','edit');
pard.Ymin.position=[4,tmin];
pard.Ymin.Width=w;
pard.Yd.object=struct('String','256','Style','edit');
pard.Yd.position=[4,td];
pard.Yd.Width=w;
pard.Ymax.object=struct('String','512','Style','edit');
pard.Ymax.position=[4,tmax];
pard.Ymax.Width=w;

pard.xyoverlapt.object=struct('String','x,y overlap (pix)','Style','text');
pard.xyoverlapt.position=[5,tp];
pard.xyoverlapt.Width=1;
pard.xyoverlap.object=struct('String','10','Style','edit');
pard.xyoverlap.position=[5,tmax];
pard.xyoverlap.Width=w;

pard.zcalc.object=struct('String','z-dependent calibration','Style','checkbox','Value',1,'Callback',{{@setvisible,obj}});
pard.zcalc.position=[6,tp];
pard.zcalc.Width=1.5;

pard.Zt.object=struct('String','Z vals (nm)','Style','text');
pard.Zt.position=[7,tp];
pard.Zt.Width=wp;
pard.Zval.object=struct('String',' 0:1000:3000','Style','edit');
pard.Zval.position=[7,tmin];
pard.Zval.Width=5-tp-wp;

pard.framerangecombinet.object=struct('String','z beyond interval (nm)','Style','text');
pard.framerangecombinet.position=[8,tp];
pard.framerangecombinet.Width=1;
pard.framerangecombine.object=struct('String','100','Style','edit');
pard.framerangecombine.position=[8,tmax];
pard.framerangecombine.Width=w;

pard.zfiltert.object=struct('String','filter beads by:','Style','text');
pard.zfiltert.position=[9,tp];
pard.zfiltert.Width=1;

pard.zfilter.object=struct('String',{{'f0','frames'}},'Style','popupmenu');
pard.zfilter.position=[9,td];
pard.zfilter.Width=2*w;


pard.fitcsplinec.object=struct('String','cspline','Style','checkbox','Callback',{{@setvisible,obj}});
pard.fitcsplinec.position=[4,1];
pard.fitcsplinec.Width=.75;

pard.fitbsplinec.object=struct('String','bspline','Style','checkbox','Callback',{{@setvisible,obj}});
pard.fitbsplinec.position=[5,1];
pard.fitbsplinec.Width=.75;

pard.roisizet.object=struct('Style','text','String','ROI (pix)'); 
pard.roisizet.position=[4,1.75];
pard.roisizet.Width=.6;
pard.roisize.object=struct('Style','edit','String','15'); 
pard.roisize.position=[4,2.35];
pard.roisize.Width=.25;

pard.roiframest.object=struct('Style','text','String','frames'); 
pard.roiframest.position=[4,2.6];
pard.roiframest.Width=.5;
pard.roiframes.object=struct('Style','edit','String','35'); 
pard.roiframes.position=[4,3.1];
pard.roiframes.Width=.25;



pard.smoothingfactort.object=struct('Style','text','String','Smoothing'); 
pard.smoothingfactort.position=[5,1.75];
pard.smoothingfactort.Width=.6;
pard.smoothingfactor.object=struct('Style','edit','String','.05'); 
pard.smoothingfactor.position=[5,2.35];
pard.smoothingfactor.Width=.25;



pard.alignzt.object=struct('Style','text','String','Align in z based on XX frames:'); 
pard.alignzt.position=[6,1.5];
pard.alignzt.Width=1.8;

pard.framewindow.object=struct('Style','edit','String','15'); 
pard.framewindow.position=[6,3.1];
pard.framewindow.Width=.25;

pard.splinepart.object=struct('String','spline smoothing','Style','text');
pard.splinepart.position=[8,2];
pard.splinepart.Width=1;
pard.splinepar.object=struct('String','0.95','Style','edit');
pard.splinepar.position=[8,3];
pard.splinepar.Width=.35;

pard.fitzc.object=struct('String','z Gauss: ','Style','checkbox','Callback',{{@setvisible,obj}});
pard.fitzc.position=[7,1];
pard.fitzc.Width=.75;
pard.fitzranget.object=struct('String','zrange nm','Style','text');
pard.fitzranget.position=[7,1.75];
pard.fitzranget.Width=.6;
pard.fitzrange.object=struct('String','-500 500','Style','edit');
pard.fitzrange.position=[7,2.35];
pard.fitzrange.Width=.5;
pard.fitzB0.object=struct('String','B0=0','Style','checkbox','Value',0);
pard.fitzB0.position=[7,2.85];
pard.fitzB0.Width=.5;


pard.fitzsxsyc.object=struct('String','sx^2-sy^2','Style','checkbox','Callback',{{@setvisible,obj}});
pard.fitzsxsyc.position=[8,1];
pard.fitzsxsyc.Width=.75;

pard.calculateZSxSy.object=struct('String','Z(Sx,Sy,X,Y,Z)','Style','checkbox','Value',0,'Callback',{{@setvisible,obj}});
pard.calculateZSxSy.position=[9,1];
pard.calculateZSxSy.Width=1.;

pard.St.object=struct('String','S (pix)','Style','text');
pard.St.position=[9,2];
pard.St.Width=wp;
pard.Smin.object=struct('String','0','Style','edit');
pard.Smin.position=[9,tmin-tp+2];
pard.Smin.Width=w;
pard.Sd.object=struct('String','.03','Style','edit');
pard.Sd.position=[9,td-tp+2];
pard.Sd.Width=w;
pard.Smax.object=struct('String','4','Style','edit');
pard.Smax.position=[9,tmax-tp+2];
pard.Smax.Width=w;

pard.inputParameters={'cam_pixelsize_nm'};
pard.plugininfo.type='ProcessorPlugin';


end
