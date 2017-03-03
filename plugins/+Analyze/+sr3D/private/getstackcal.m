function [SXY,beadgood]=getstackcal(beadsh,p,X,Y,axall)
global stackcal_testfit

alignzf0=p.alignz;
alignxcorr=p.alignzxcorr;

zc=p.spatialcalibration & p.zcalc &p.beaddistribution.Value==2; %z-dependent calibration

sstack=size(beadsh(1).stack.image);

%for sorting

fmedian=nanmedian([beadsh(:).f0]);
halfstoreframes=(size(beadsh(1).stack.image,3)-1)/2;
if zc||alignzf0||alignxcorr
    df=0;
    fzero=halfstoreframes+1; %only used for plotting. Take out? f0 should be in center now.
else
    f0=[beadsh(:).f0];
    fzero=round(fmedian);
    df=(f0-fmedian);
end
 psfx=[beadsh(:).psfx0];psfy=[beadsh(:).psfy0];
dpsfx=(psfx-median(psfx(~isnan(psfx))))*10;
dpsfy=(psfy-median(psfy(~isnan(psfy))))*10;
 
for Z=1:length(p.Zrange)-1
    allstacks=zeros(sstack(1),sstack(2),sstack(3),length(beadsh))+NaN;
    for B=length(beadsh):-1:1
        beadz0=(beadsh(B).f0-p.fglass(beadsh(B).filenumber))*p.dz;
%         
%         if p.alignz||p.beaddistribution.Value==2
%              beadz=(beadsh(B).loc.frame*p.dz)-beadz0;
%         else 
%             beadz=(beadsh(B).loc.frame)*p.dz;
%         end
%         zglass=p.fglass*p.dz;
       
        if zc
            if p.zfilter.Value==1 % f0

                if beadz0>p.Zrange(Z)-p.framerangecombine*p.dz && beadz0<p.Zrange(Z+1)+p.framerangecombine
%                    sx=beadsh(B).loc.PSFxpix; 
                   allstacks(:,:,:,B)=beadsh(B).stack.image;
                end             
            else
                zs=(beadsh(B).stack.framerange-p.fglass(beadsh(B).filenumber))*p.dz;
                goodz=zs>p.Zrange(Z)-p.framerangecombine & zs<p.Zrange(Z+1)+p.framerangecombine;
                allstacks(:,:,goodz,B)=beadsh(B).stack.image(:,:,goodz);
            end
        else
            allstacks(:,:,:,B)=beadsh(B).stack.image;
            
        end
        
        stackh=allstacks(:,:,:,B);
%         inth=nansum(stackh(:)); %take into account when sorting
        goodvs(B)=sum(~isnan(stackh(:)))/numel(stackh);
   
        if  halfstoreframes<length(beadsh(B).stack.framerange)
            dframe(B)=beadsh(B).stack.framerange(halfstoreframes+1)-beadsh(B).f0;
        else
            dframe(B)=0;
        end
    end
    
    %sort
    devs=(df.^2+dpsfx.^2+dpsfy.^2)./goodvs;
    
    if p.alignzxcorr
        fw2=round((p.framewindow-1)/2);
    else
        fw2=2;
    end
    
    tgt=[num2str(X) num2str(Y) num2str(Z)];
    ax=maketgax(axall.hspline_scatter,tgt);  
    ax2=maketgax(axall.hspline_psf,tgt);
    
 %get calibrations
    SXY(Z)=getsxyinit(p,X,Y,Z);
    
    [~,sortinddev]=sort(devs);
    allrois=allstacks(:,:,:,sortinddev);
    
    if ~alignxcorr&&alignzf0
        zshift=dframe(sortinddev)-median(dframe);
    else
        zshift=[];
    end
%     zshift=-zshift;
    midrange=halfstoreframes+round(median(dframe));
    
    [corrPSF,shiftedstack,shift,beadgood]=registerPSF3D(allrois,struct('framerange',midrange+1-fw2:midrange+1+fw2,'alignz',alignxcorr,'zshiftf0',zshift),{ax, ax2});

        %undo sorting by deviation to associate beads again to their
    %bead number
  [~,sortback]=sort(sortinddev);
    shiftedstack=shiftedstack(:,:,:,sortback);
    shift=shift(sortback,:);
    beadgood=beadgood(sortback);
    allrois=allstacks;
  
     %calculate dual-color overlay for display
            allroisol=zeros(size(allrois,1),size(allrois,2),size(allrois,3),size(allrois,4),3);
            for k=1:size(allrois,4)
                ssh=shiftedstack(:,:,:,k);
                allroisol(:,:,:,k,1)=ssh/nansum(ssh(:));
                allroisol(:,:,:,k,2)=corrPSF/nansum(corrPSF(:));
%                 allroisol(:,:,:,k,3)=0.5*(allroisol(:,:,:,k,1)+allroisol(:,:,:,k,2));
            end
            axallps=maketgax(axall.hspline_overlay,tgt);
%             axallps=obj.initaxis('overlayPSF');
            imageslicer(allroisol,'rgb',1,'Parent',axallps.Parent)
%             if p.alignz
%                 ax=obj.initaxis('zshift');scatter3(beadpos(:,1),beadpos(:,2),shift(:,3))
%             end
%             drawnow
            
            %cut out the central part of the PSF correspoinding to the set
            %Roisize in x,y and z
            
            [im,ind]=nanmax(corrPSF(:));
            
            if ind>1
            [x,y,z]=ind2sub(size(corrPSF),ind);
            dRx=round((p.roisize-1)/2);
            dz=round((p.roiframes-1)/2);
            rangex=x-dRx:x+dRx;
            rangey=y-dRx:y+dRx;
            rangez=max(1,z-dz):min(size(corrPSF,3),z+dz);
            
            %normalize PSF
            
            minPSF=min(corrPSF(:));
            corrPSFn=corrPSF-minPSF;
            intglobal=nanmean(nansum(nansum(corrPSFn(rangex,rangey,z-1:z+1),1),2));
            for k=1:size(corrPSF,3)
                corrPSFn(:,:,k)=corrPSFn(:,:,k)/intglobal;
                
%                 int(k)=nansum(nansum(corrPSFn(:,:,k)));
%                 if int(k)>0
% %                 corrPSFn(:,:,k)=corrPSFn(:,:,k)/int(k);
%                 end
            end
            corrPSFn(isnan(corrPSFn))=0;
            
            corrPSFs=corrPSFn(rangex,rangey,rangez);
            PSFgood=true;
            else
                corrPSFn=0*corrPSF;
                corrPSFs=zeros(p.roisize,p.roisize,p.roiframes);
                PSFgood=false;
            end
            
            
            %Normalization of frames??? where to do that? berfore b-spline?
            %calculatae b-spline coefficients
            
            %calculate effective smoothing factor
            if p.smooth
            lambda=p.smoothingfactor*100;
            if length(lambda)<2
                lambda(2:3)=lambda;
            elseif length(lambda)<3
                lambda=[lambda(1) lambda(1) lambda(2)];
            end
            lambda(1:2)=lambda(1:2)/p.cam_pixelsize_nm;
            lambda(3)=lambda(3)/p.dz;
            else
                lambda=0;
            end
            
            b3_0=bsarray(double(corrPSFs),'lambda',lambda);
            
            %calculate double-sampled PSF for c-spline
            
            zhd=1:1:b3_0.dataSize(3);
%             if p.doublesample
%                 dxxhd=0.5;
 dxxhd=1;
                [XX,YY,ZZ]=meshgrid(1:dxxhd:b3_0.dataSize(1),1:dxxhd:b3_0.dataSize(2),zhd);
                corrPSFhd = interp3_0(b3_0,XX,YY,ZZ,0);
%             else
               
%                 corrPSFhd=corrPSFs;
%             end
                
  
            
            if p.fitcsplinec
            tic
            spline = Spline3D_v2(corrPSFhd);
            coeff = spline.coeff;
%             [np_psf,coeff]=generate_psf_to_spline_YLJ(corrPSFhd,p.roisize,dz+1);
            toc
%             bspline
%             np_psf = np_psf/max(np_psf(:));
%             for i = 1:size(np_psf,3)
%                 np_psf(:,:,i) = np_psf(:,:,i)/sum(sum(np_psf(:,:,i)));
%             end
            else
                coeff=0;
            end

            
%             save
            bspline.bslpine=b3_0;
            bspline.isEM=p.EMon;
            cspline.coeff=coeff;
            cspline.isEM=p.EMon;
%             cspline.doublesample=p.doublesample;
            cspline.z0=round((b3_0.dataSize(3)+1)/2);
            cspline.dz=p.dz;
            bspline.z0=round((b3_0.dataSize(3)+1)/2);
            bspline.dz=p.dz;            
%             save([path filesep file],'bspline','cspline','z0','dz')
%     if p.fitbsplinec      
        SXY(Z).splinefit.bspline=bspline;
        SXY(Z).splinefit.PSF=corrPSF;
         SXY(Z).splinefit.PSFsmooth=corrPSFhd;
%     end
    if p.fitcsplinec
        SXY(Z).splinefit.cspline=cspline;
    end
    
    if PSFgood
            %plot graphs
            ax=maketgax(axall.hspline_psfz,tgt);
%             ax=obj.initaxis('PSFz');
             framerange0=max(p.fminmax(1),fzero-halfstoreframes):min(fzero+halfstoreframes,p.fminmax(2));
             halfroisizebig=(size(shiftedstack,1)-1)/2;
            ftest=fzero-framerange0(1)+1;
            
            xt=halfroisizebig+1;
            yt=halfroisizebig+1;
            
            ftest=z;
            xt=x;
            yt=y;
            zpall=squeeze(shiftedstack(xt,yt,:,beadgood));
            zpall2=squeeze(allrois(xt,yt,:,beadgood));
            xpall=squeeze(shiftedstack(:,yt,ftest,beadgood));
            xpall2=squeeze(allrois(:,yt,ftest,beadgood));
            
            for k=1:size(zpall,2)
                zpall(:,k)=zpall(:,k)/nanmax(zpall(:,k));
                xpall(:,k)=xpall(:,k)/nanmax(xpall(:,k));
                zpall2(:,k)=zpall2(:,k)/nanmax(zpall2(:,k));
                xpall2(:,k)=xpall2(:,k)/nanmax(xpall2(:,k));                
            end
            
            zprofile=squeeze(corrPSFn(xt,yt,:));
            zprofile=zprofile/max(zprofile);
            mphd=round((size(corrPSFhd,1)+1)/2);
            zprofilehd=squeeze(corrPSFhd(mphd,mphd,:));
            zprofilehd=zprofilehd/max(zprofilehd);            
            
            xprofile=squeeze(corrPSFn(:,yt,ftest));
            xprofile=xprofile/max(xprofile);
            mpzhd=round((size(corrPSFhd,3)+1)/2+1);
            xprofilehd=squeeze(corrPSFhd(:,mphd,mpzhd));
            xprofilehd=xprofilehd/max(xprofilehd);
            
            dxxx=0.1;
            xxx=1:dxxx:b3_0.dataSize(1);
            zzzt=0*xxx+mpzhd;
%             zzzt=0*xxx+ftest-rangez(1)-0*framerange0(1)+1;
%             xbs= interp3_0(b3_0,xxx,0*xxx+b3_0.dataSize(1)/2+.5,zzzt);
            xbs= interp3_0(b3_0,0*xxx+b3_0.dataSize(1)/2+.5,xxx,zzzt);
            xbs=xbs/max(xbs);
            
            zzz=1:dxxx:b3_0.dataSize(3);xxxt=0*zzz+b3_0.dataSize(1)/2+.5;
            zbs= interp3_0(b3_0,xxxt,xxxt,zzz); 
            zbs=zbs/max(zbs);

           hold(ax,'off')
             plot(ax,framerange0,zpall2(1:length(framerange0),:),'m:')
            hold(ax,'on')
             plot(ax,framerange0,zpall(1:length(framerange0),:),'c')
            plot(ax,framerange0',zprofile(1:length(framerange0)),'k*')
            plot(ax,rangez,zprofilehd,'ko')
            plot(ax,zzz+rangez(1)+framerange0(1)-2,zbs,'b','LineWidth',2)
            
            xrange=-halfroisizebig:halfroisizebig;  
            
            
            ax=maketgax(axall.hspline_psfx,tgt);
%              ax=obj.initaxis('PSFx');
            hold(ax,'off')
            plot(ax,xrange,xpall2,'m:')
            hold(ax,'on')
            plot(ax,xrange,xpall,'c')
            plot(ax,xrange,xprofile,'k*-')
            plot(ax,xrange(rangex),xprofilehd,'ko')
            plot(ax,(xxx-(b3_0.dataSize(1)+1)/2),xbs,'b','LineWidth',2)
            
            drawnow
            %quality control: refit all beads
            ax=maketgax(axall.hspline_validate,tgt);
%             hold(ax,'off')
%             axes(ax)
            drawnow
            if isempty(stackcal_testfit)||stackcal_testfit
                testallrois=allrois(:,:,:,beadgood);
                testallrois(isnan(testallrois))=0;
            testfit(testallrois,cspline.coeff,p,{},ax)
            testfit(corrPSF,cspline.coeff,p,{'k','LineWidth',2},ax)
            end
   end 


end
end


function testfit(teststack,coeff,p,linepar,ax)
if nargin<4
    linepar={};
elseif ~iscell(linepar)
    linepar={linepar};
end
d=round((size(teststack,1)-p.roisize)/2);
            range=d+1:d+p.roisize;
            
% if p.doublesample
%     coefffak=4;
%     fitterGPU=@GPUmleFit_LM_v2;
%     fitterCPU=@kernel_MLEfit_Spline_LM_SMAP_v2;
% else
    coefffak=1;
    fitterGPU=@callYimingFitter;
    fitterCPU=@kernel_MLEfit_Spline_LM_SMAP_v2_nointerp;
% end  
    for k=1:size(teststack,4)
        try
          [P] =  fitterGPU(single(squeeze(teststack(range,range,:,k))),single(coefffak*coeff),100,5,0);
        catch err
%             err
            disp('run on CPU')        
            [P] =  fitterCPU(teststack(range,range,:,k),(coefffak*coeff),single(p.roisize),100);
        end
    z=(1:size(P,1))'-1;
    plot(ax,z,P(:,5),linepar{:})
    hold(ax,'on')
    end
end