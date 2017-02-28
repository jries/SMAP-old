function [SXY,beadgood]=getstackcal(beadsh,p,X,Y,axall)
global stackcal_testfit
zc=p.spatialcalibration & p.zcalc &p.beaddistribution.Value==2;

sstack=size(beadsh(1).stack.image);

%for sorting

fmedian=nanmedian([beadsh(:).f0]);
halfstoreframes=(size(beadsh(1).stack.image,3)-1)/2;
if zc
    df=0;
    fzero=halfstoreframes+1;
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
        beadz0=(beadsh(B).f0-p.fglass)*p.dz;
%         
%         if p.alignz||p.beaddistribution.Value==2
%              beadz=(beadsh(B).loc.frame*p.dz)-beadz0;
%         else 
%             beadz=(beadsh(B).loc.frame)*p.dz;
%         end
%         zglass=p.fglass*p.dz;
       
        if zc
            if p.zfilter.Value==1 % f0

                if beadz0>p.Zrange(Z)-p.framerangecombine && beadz0<p.Zrange(Z+1)+p.framerangecombine
                   sx=beadsh(B).loc.PSFxpix; 
                   allstacks(:,:,:,B)=beadsh(B).stack.image;
                end             
            else
                zs=(beadsh(B).stack.framerange-p.fglass)*p.dz;
                goodz=zs>p.Zrange(Z)-p.framerangecombine & zs<p.Zrange(Z+1)+p.framerangecombine;
                allstacks(:,:,goodz,B)=beadsh(B).stack.image(:,:,goodz);
            end
        else
            allstacks(:,:,:,B)=beadsh(B).stack.image;
            
        end
        
        stackh=allstacks(:,:,:,B);
%         inth=nansum(stackh(:)); %take into account when sorting
        goodvs(B)=sum(~isnan(stackh(:)))/numel(stackh);
        
    end
    
    %sort
    devs=(df.^2+dpsfx.^2+dpsfy.^2)./goodvs;
    
    fw2=round((p.framewindow-1)/2);
    
    tgt=[num2str(X) num2str(Y) num2str(Z)];
    ax=maketgax(axall.hspline_scatter,tgt);  
    ax2=maketgax(axall.hspline_psf,tgt);
    
 %get calibrations
    SXY(Z)=getsxyinit(p,X,Y,Z);
    
    [~,sortinddev]=sort(devs);
    allrois=allstacks(:,:,:,sortinddev);
    
    [corrPSF,shiftedstack,shift,beadgood]=registerPSF3D(allrois,struct('framerange',halfstoreframes+1-fw2:halfstoreframes+1+fw2,'alignz',zc|p.alignz),{ax, ax2});

        %undo sorting by deviation to associate beads again to their
    %bead number
    numbers=1:length(sortinddev);
    shiftedstack=shiftedstack(:,:,:,numbers(sortinddev));
    shift=shift(numbers(sortinddev),:);
    beadgood=beadgood(numbers(sortinddev));

    
     %calculate dual-color overlay for display
            allroisol=zeros(size(allrois,1),size(allrois,2),size(allrois,3),size(allrois,4),3);
            for k=1:size(allrois,4)
                ssh=shiftedstack(:,:,:,k);
                allroisol(:,:,:,k,1)=ssh/nansum(ssh(:));
                allroisol(:,:,:,k,2)=corrPSF/nansum(corrPSF(:));
                allroisol(:,:,:,k,3)=0.5*(allroisol(:,:,:,k,1)+allroisol(:,:,:,k,2));
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
            if im>1
            [x,y,z]=ind2sub(size(corrPSF),ind);
            dRx=round((p.roisize-1)/2);
            dz=round((p.roiframes-1)/2);
            rangex=x-dRx:x+dRx;
            rangey=y-dRx:y+dRx;
            rangez=z-dz:z+dz;
            
            %normalize PSF
            
            minPSF=min(corrPSF(:));
            corrPSFn=corrPSF-minPSF;
            corrPSFn(isnan(corrPSFn))=0;
            for k=1:size(corrPSF,3)
                int=nansum(nansum(corrPSFn(:,:,k)));
                if int>0
                corrPSFn(:,:,k)=corrPSFn(:,:,k)/int;
                end
            end
            
            
            corrPSFs=corrPSFn(rangex,rangey,rangez);
            PSFgood=true;
            else
                corrPSFn=0*corrPSF;
                corrPSFs=zeros(p.roisize,p.roisize,p.roiframes);
                PSFgood=false;
            end
            
            
            %Normalization of frames??? where to do that? berfore b-spline?
            %calculatae b-spline coefficients
            b3_0=bsarray(double(corrPSFs),'lambda',p.smoothingfactor);
            
            %calculate double-sampled PSF for c-spline
            
            zhd=1:1:b3_0.dataSize(3);
%             if p.doublesample
%                 dxxhd=0.5;
%                 [XX,YY,ZZ]=meshgrid(1:dxxhd:b3_0.dataSize(1),1:dxxhd:b3_0.dataSize(2),zhd);
%                 corrPSFhd = interp3_0(b3_0,XX,YY,ZZ,0);
%             else
                dxxhd=1;
                corrPSFhd=corrPSFs;
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
    if p.fitbsplinec      
        SXY(Z).splinefit.bspline=bspline;
    end
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
            zpall=squeeze(shiftedstack(halfroisizebig+1,halfroisizebig+1,:,beadgood));
            zpall2=squeeze(allrois(halfroisizebig+1,halfroisizebig+1,:,beadgood));
            xpall=squeeze(shiftedstack(:,halfroisizebig+1,ftest,beadgood));
            xpall2=squeeze(allrois(:,halfroisizebig+1,ftest,beadgood));
            
            for k=1:size(zpall,2)
                zpall(:,k)=zpall(:,k)/max(zpall(:,k));
                xpall(:,k)=xpall(:,k)/max(xpall(:,k));
                zpall2(:,k)=zpall2(:,k)/max(zpall2(:,k));
                xpall2(:,k)=xpall2(:,k)/max(xpall2(:,k));                
            end
            
            zprofile=squeeze(corrPSFn(halfroisizebig+1,halfroisizebig+1,:));
            zprofile=zprofile/max(zprofile);
            mphd=round((size(corrPSFhd,1)+1)/2);
            zprofilehd=squeeze(corrPSFhd(mphd,mphd,:));
            zprofilehd=zprofilehd/max(zprofilehd);            
            
            xprofile=squeeze(corrPSFn(:,halfroisizebig+1,ftest));
            xprofile=xprofile/max(xprofile);
            mpzhd=round(size(corrPSFhd,3)+1)/2+1;
            xprofilehd=squeeze(corrPSFhd(:,mphd,mpzhd));
            xprofilehd=xprofilehd/max(xprofilehd);
            
            dxxx=0.1;
            xxx=1:dxxx:b3_0.dataSize(1);zzzt=0*xxx+ftest-rangez(1)-0*framerange0(1)+1;
%             xbs= interp3_0(b3_0,xxx,0*xxx+b3_0.dataSize(1)/2+.5,zzzt);
            xbs= interp3_0(b3_0,0*xxx+b3_0.dataSize(1)/2+.5,xxx,zzzt);
            xbs=xbs/max(xbs);
            
            zzz=1:dxxx:b3_0.dataSize(3);xxxt=0*zzz+b3_0.dataSize(1)/2+.5;
            zbs= interp3_0(b3_0,xxxt,xxxt,zzz); 
            zbs=zbs/max(zbs);

            hold off
             plot(framerange0,zpall2(1:length(framerange0)),'m:')
             hold on
             plot(framerange0,zpall(1:length(framerange0)),'c')
            plot(framerange0',zprofile(1:length(framerange0)),'k*')
            plot(zhd+rangez(1)+framerange0(1)-2,zprofilehd,'ko')
            plot(zzz+rangez(1)+framerange0(1)-2,zbs,'b','LineWidth',2)
            
            xrange=-halfroisizebig:halfroisizebig;  
            
            
            ax=maketgax(axall.hspline_psfx,tgt);
%              ax=obj.initaxis('PSFx');
            hold off
            plot(xrange,xpall2,'m:')
            hold on
            plot(xrange,xpall,'c')
            plot(xrange,xprofile,'k*-')
            plot((-mphd+1:mphd-1)/dxxhd,xprofilehd,'ko')
            plot((xxx-(b3_0.dataSize(1)+1)/2),xbs,'b','LineWidth',2)
            
            drawnow
            %quality control: refit all beads
            ax=maketgax(axall.hspline_validate,tgt);
            hold off
    
            if isempty(stackcal_testfit)||stackcal_testfit
            testfit(allrois(:,:,:,beadgood),cspline.coeff,p)
            testfit(corrPSF,cspline.coeff,p,{'k','LineWidth',2})
            end
   end 


end
end


function testfit(teststack,coeff,p,linepar)
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
    fitterGPU=@GPUmleFit_LM_noInterp;
    fitterCPU=@kernel_MLEfit_Spline_LM_SMAP_v2_nointerp;
% end  
    for k=1:size(teststack,4)
        try
          [P] =  fitterGPU(single(squeeze(teststack(range,range,:,k))),single(coefffak*coeff),100,5,0);
        catch err
%             err
            disp('run on CPU')        
            [P] =  fitterCPU(teststack(range,range,:,k),(4*coeff),single(p.roisize),100);
        end
    z=(1:size(P,1))'-1;
    plot(z,P(:,5),linepar{:})
    hold on
    end
end