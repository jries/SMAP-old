function [splinefit,indgood]=getstackcal_so(beads,p)
global stackcal_testfit
isastig=contains(p.modality,'astig');
alignzastig=isastig&contains(p.modality,'astig');
zcorr=contains(p.zcorr,'corr');

% iter=p.iter;
% alignzf0=p.alignz;
% alignxcorr=p.alignzxcorr;

% zc=p.spatialcalibration & p.zcalc &p.beaddistribution.Value==2; %z-dependent calibration

sstack=size(beads(1).stack.image);

%for sorting
% for Z=1:length(p.Zrange)-1
%     beadshz=beadsh(indgood{Z});
    
%     fmedian=nanmedian([beads(:).f0]);
    halfstoreframes=(size(beads(1).stack.image,3)-1)/2;
%     if zc||alignzf0||alignxcorr
%         df=0;
%         fzero=halfstoreframes+1; %only used for plotting. Take out? f0 should be in center now.
%     else
%         f0=[beadshz(:).f0];
%         fzero=round(fmedian);
%         df=(f0-fmedian);
%     end
    if isastig    
        for B=length(beads):-1:1
            if  halfstoreframes<length(beads(B).stack.framerange)
                dframe(B)=beads(B).stack.framerange(halfstoreframes+1)-beads(B).f0;
            else
                dframe(B)=NaN;
            end
        end
            %remove outliers:
        badind=abs(dframe-nanmedian(dframe))>10|isnan(dframe);
    %     beads(badind)=[];
    %     dframe(badind)=[];
    %     goodvs(badind)=[];
        beads(badind)=[];
    

        psfx=[beads(:).psfx0];psfy=[beads(:).psfy0];
        dpsfx=(psfx-median(psfx(~isnan(psfx))))*10;
        dpsfy=(psfy-median(psfy(~isnan(psfy))))*10;
    else
        dpsfx=0;dpsfy=0;
    end
    

    allstacks=zeros(sstack(1),sstack(2),sstack(3),length(beads))+NaN;
    goodvs=[];
    for B=length(beads):-1:1
        allstacks(:,:,:,B)=beads(B).stack.image;
        stackh=allstacks(:,:,:,B);
        goodvs(B)=sum(~isnan(stackh(:)))/numel(stackh);
%         if  halfstoreframes<length(beads(B).stack.framerange)
%             dframe(B)=beads(B).stack.framerange(halfstoreframes+1)-beads(B).f0;
%         else
%             dframe(B)=0;
%         end
    end
    

    
    mstack=nanmean(allstacks,4);
    mstack=mstack-nanmin(mstack(:));
    mstack=mstack/nansum(mstack(:));
    for k=length(beads):-1:1
    	stackh=(allstacks(:,:,:,k));
        stackh=stackh-nanmin(stackh(:));
        stackh=stackh/nansum(stackh(:));
        dstack(k)=sum((stackh(:)-mstack(:)).^2);
    end
    dstack=dstack/mean(dstack);    
    devs=(dpsfx.^2+dpsfy.^2+dstack)./goodvs;

    if zcorr
        fw2=round((p.zcorrframes-1)/2);
    else
        fw2=2;
    end
    
    ax=axes('Parent',uitab(p.tabgroup,'Title','scatter'));
%     ax2=axes('Parent',uitab(p.tabgroup,'Title','PSF'));
%     ax=maketgax(axall.hspline_scatter,tgt);  
%     ax2=maketgax(axall.hspline_psf,tgt);
    

    [~,sortinddev]=sort(devs);
    allrois=allstacks(:,:,:,sortinddev);
    
    if alignzastig
        zshift=dframe(sortinddev)-round(median(dframe));
    else
        zshift=[];
    end
    zshift=-zshift;
    midrange=halfstoreframes;%+round(median(dframe))+1;
    
    [corrPSF,shiftedstack,shift,beadgood]=registerPSF3D_so(allrois,struct('framerange',midrange-fw2:midrange+fw2,'alignz',zcorr,'zshiftf0',zshift,'beadfilterf0',false),{ax});

        %undo sorting by deviation to associate beads again to their
    %bead number
    [~,sortback]=sort(sortinddev);
    shiftedstack=shiftedstack(:,:,:,sortback);
    shift=shift(sortback,:);
    beadgood=beadgood(sortback);

    indgood=beadgood;
    allrois=allstacks;
  
     %calculate dual-color overlay for display
        allroisol=zeros(size(allrois,1),size(allrois,2),size(allrois,3),size(allrois,4),3);
        psfnorm=corrPSF/nansum(corrPSF(:));
        for k=1:size(allrois,4)
            ssh=shiftedstack(:,:,:,k);
            allroisol(:,:,:,k,1)=ssh/nansum(ssh(:));
            allroisol(:,:,:,k,2)=psfnorm;
%                 allroisol(:,:,:,k,3)=0.5*(allroisol(:,:,:,k,1)+allroisol(:,:,:,k,2));
        end
%         axallps=maketgax(axall.hspline_overlay,tgt);
%             axallps=obj.initaxis('overlayPSF');
%         imageslicer(allroisol,'rgb',1,'Parent',axallps.Parent)
%             if p.alignz
%                 ax=obj.initaxis('zshift');scatter3(beadpos(:,1),beadpos(:,2),shift(:,3))
%             end
%             drawnow

        %cut out the central part of the PSF correspoinding to the set
        %Roisize in x,y and z

        scorrPSF=size(corrPSF);
        x=round((scorrPSF(1)+1)/2);y=round((scorrPSF(2)+1)/2);
%         if strcmp(p.modality,'astigmatic') %align to brightest 
%             [im,ind]=nanmax(corrPSF(:));  
%             if ind>1
%                 [x,y,z]=ind2sub(size(corrPSF),ind);
%             end
%         end

        dRx=round((p.ROIxy-1)/2);
        dz=round((p.ROIz-1)/2);
        rangex=x-dRx:x+dRx;
        rangey=y-dRx:y+dRx;

        z=midrange;%always same reference: z=f0


        rangez=max(1,z-dz):min(size(corrPSF,3),z+dz);
%             framess=
        %normalize PSF
        centpsf=corrPSF(2:end-1,2:end-1,2:end-1); %cut out rim from shift
        minPSF=min(centpsf(:));
        corrPSFn=corrPSF-minPSF;
        intglobal=nanmean(nansum(nansum(corrPSFn(rangex,rangey,z-1:z+1),1),2));
        for k=1:size(corrPSF,3)
            corrPSFn(:,:,k)=corrPSFn(:,:,k)/intglobal;

%                 int(k)=nansum(nansum(corrPSFn(:,:,k)));
%                 if int(k)>0
% %                 corrPSFn(:,:,k)=corrPSFn(:,:,k)/int(k);
%                 end
        end
        shiftedstack=shiftedstack/intglobal;
        corrPSFn(isnan(corrPSFn))=0;
        corrPSFn(corrPSFn<0)=0;
        corrPSFs=corrPSFn(rangex,rangey,rangez);
        
        PSFgood=true;
%             else
%                 corrPSFn=0*corrPSF;
%                 corrPSFs=zeros(p.roisize,p.roisize,p.roiframes);
%                 PSFgood=false;
%             end
            
            
        %Normalization of frames??? where to do that? berfore b-spline?
        %calculatae b-spline coefficients

        %calculate effective smoothing factor
        
        lambdax=p.smoothxy/p.cam_pixelsize_um(1)/100;
        lambdaz=p.smoothz/p.dz;
%         if length(lambda)<2
%             lambda(2:3)=lambda;
%         elseif length(lambda)<3
%             lambda=[lambda(1) lambda(1) lambda(2)];
%         end
%         lambda(1:2)=lambda(1:2).;
%         lambda(3)=lambda(3)/p.dz;
%         else
%             lambda=0;
%         end
        lambda=[lambdax lambdax lambdaz];
        b3_0=bsarray(double(corrPSFs),'lambda',lambda);

        %calculate double-sampled PSF for c-spline

        zhd=1:1:b3_0.dataSize(3);
%             if p.doublesample
%                 dxxhd=0.5;
        dxxhd=1;
        [XX,YY,ZZ]=meshgrid(1:dxxhd:b3_0.dataSize(1),1:dxxhd:b3_0.dataSize(2),zhd);
        corrPSFhd = interp3_0(b3_0,XX,YY,ZZ,0);
        
%         psfshow(:,:,:,1)=corrPSFs;
%         psfshow(:,:,:,2)=corrPSFhd;
%         psfshow(1,1,1,3)=0;
%         imageslicer(psfshow,'Parent',ax2,'Title','PSF: 1.average, 2.smooth')
%             else
               
%                 corrPSFhd=corrPSFs;
%             end
                
  

%         if p.fitcsplinec
        tic
%             profile on
        spline = Spline3D_v2(corrPSFhd);
%             profile viewer
        coeff = spline.coeff;
%             [np_psf,coeff]=generate_psf_to_spline_YLJ(corrPSFhd,p.roisize,dz+1);
        toc
%             bspline
%             np_psf = np_psf/max(np_psf(:));
%             for i = 1:size(np_psf,3)
%                 np_psf(:,:,i) = np_psf(:,:,i)/sum(sum(np_psf(:,:,i)));
%             end
%         else
%             coeff=0;
%         end

            
%             save
        bspline.bslpine=b3_0;
%         bspline.isEM=p.EMon;
        cspline.coeff=coeff;
%         cspline.isEM=p.EMon;
%             cspline.doublesample=p.doublesample;
        cspline.z0=round((b3_0.dataSize(3)+1)/2);
        cspline.dz=p.dz;
        bspline.z0=round((b3_0.dataSize(3)+1)/2);
        bspline.dz=p.dz;            
%             save([path filesep file],'bspline','cspline','z0','dz')
%     if p.fitbsplinec      
        splinefit.bspline=bspline;
        splinefit.PSF=corrPSF;
        splinefit.PSFsmooth=corrPSFhd;
%     end
 
        splinefit.cspline=cspline;

    
        if PSFgood
            %plot graphs
            ax=axes(uitab(p.tabgroup,'Title','PSFz'));
%             ax=obj.initaxis('PSFz');
             framerange0=max(p.fminmax(1)):min(p.fminmax(2));
             halfroisizebig=(size(shiftedstack,1)-1)/2;
%             ftest=fzero-framerange0(1)+1;
            
%             xt=halfroisizebig+1;
%             yt=halfroisizebig+1;
            
            ftest=z;
            xt=x;
            yt=y;
            zpall=squeeze(shiftedstack(xt,yt,:,beadgood));
            zpall2=squeeze(allrois(xt,yt,:,beadgood));
            xpall=squeeze(shiftedstack(:,yt,ftest,beadgood));
            xpall2=squeeze(allrois(:,yt,ftest,beadgood));
            
            for k=1:size(zpall,2)
%                 zpall(:,k)=zpall(:,k)/nanmax(zpall(:,k));
%                 xpall(:,k)=xpall(:,k)/nanmax(xpall(:,k));
                zpall2(:,k)=zpall2(:,k)/nanmax(zpall2(:,k));
                xpall2(:,k)=xpall2(:,k)/nanmax(xpall2(:,k));                
            end
            
            zprofile=squeeze(corrPSFn(xt,yt,:));
%             zprofile=zprofile/max(zprofile);
            mphd=round((size(corrPSFhd,1)+1)/2);
            zprofilehd=squeeze(corrPSFhd(mphd,mphd,:));
%             zprofilehd=zprofilehd/max(zprofilehd);            
            
            xprofile=squeeze(corrPSFn(:,yt,ftest));
%             xprofile=xprofile/max(xprofile);
            mpzhd=round((size(corrPSFhd,3)+1)/2+1);
            xprofilehd=squeeze(corrPSFhd(:,mphd,mpzhd));
%             xprofilehd=xprofilehd/max(xprofilehd);
            
            dxxx=0.1;
            xxx=1:dxxx:b3_0.dataSize(1);
            zzzt=0*xxx+mpzhd;
%             zzzt=0*xxx+ftest-rangez(1)-0*framerange0(1)+1;
%             xbs= interp3_0(b3_0,xxx,0*xxx+b3_0.dataSize(1)/2+.5,zzzt);
            xbs= interp3_0(b3_0,0*xxx+b3_0.dataSize(1)/2+.5,xxx,zzzt);
%             xbs=xbs/max(xbs);
            
            zzz=1:dxxx:b3_0.dataSize(3);xxxt=0*zzz+b3_0.dataSize(1)/2+.5;
            zbs= interp3_0(b3_0,xxxt,xxxt,zzz); 
%             zbs=zbs/max(zbs);

           hold(ax,'off')
%              plot(ax,framerange0,zpall2(1:length(framerange0),:),'m:')
            
             plot(ax,framerange0,zpall(1:length(framerange0),:),'c')
             hold(ax,'on')
            plot(ax,framerange0',zprofile(1:length(framerange0)),'k*')
            plot(ax,rangez+framerange0(1)-1,zprofilehd,'ko')
            plot(ax,zzz+rangez(1)+framerange0(1)-2,zbs,'b','LineWidth',2)
            
            xrange=-halfroisizebig:halfroisizebig;  
            
             ax=axes(uitab(p.tabgroup,'Title','PSFx'));
%             ax=maketgax(axall.hspline_psfx,tgt);
%              ax=obj.initaxis('PSFx');
            hold(ax,'off')
%             plot(ax,xrange,xpall2,'m:')
            
            plot(ax,xrange,xpall,'c')
            hold(ax,'on')
            plot(ax,xrange,xprofile,'k*-')
            plot(ax,xrange(rangex),xprofilehd,'ko')
            plot(ax,(xxx-(b3_0.dataSize(1)+1)/2),xbs,'b','LineWidth',2)
            
            drawnow
            %quality control: refit all beads
             ax=axes(uitab(p.tabgroup,'Title','validate'));
%             ax=maketgax(axall.hspline_validate,tgt);
%             hold(ax,'off')
%             axes(ax)
            drawnow
            if isempty(stackcal_testfit)||stackcal_testfit
                testallrois=allrois(:,:,:,beadgood);
                testallrois(isnan(testallrois))=0;
                zall=testfit(testallrois,cspline.coeff,p,{},ax);
                corrPSFfit=corrPSF/max(corrPSF(:))*max(testallrois(:)); %bring back to some reasonable photon numbers;
                zref=testfit(corrPSFfit,cspline.coeff,p,{'k','LineWidth',2},ax);

                drawnow
                 ax=axes(uitab(p.tabgroup,'Title','stripes'));
                teststripes(cspline.coeff,p,ax);
            end
%             ztest=zall;
%             for k=1:size(ztest,2)
%                 xa=find(ztest(:,k)>15.5,1,'first');
%                 xa=19;
%                 xb=xa+1;
%                 ya=ztest(xa,k);
%                 yb=ztest(xb,k);
%                 y0=16;
%                 framezero(k)=(xb-xa)/(yb-ya)*(y0-ya)+xa;
%             end
%                 xa=find(zref>15.5,1,'first');
%                 xa=19;
%                 xb=xa+1;
%                 ya=zref(xa);
%                 yb=zref(xb);
%                 y0=16;
%                 framezeropsf=(xb-xa)/(yb-ya)*(y0-ya)+xa;
        end 


% end
end

function teststripes(coeff,p,ax)
tt=tic;

zr=0:0.2:p.ROIz;
xr=0:0.05:p.ROIxy;
hz=zeros(1,length(zr)-1);
hx=zeros(1,length(xr)-1);
hy=hx;
while toc(tt)<30
    nn=rand(11,11,10000,'single');
    P=callYimingFitter(nn,single(coeff),50,5,0,0);
    
    hz=histcounts(P(:,5),zr)+hz;
    hx=histcounts(P(:,1),xr)+hx;
    hy=histcounts(P(:,2),xr)+hy;
    
end

hz(1)=[];hz(end)=[];
hz(1)=0;hz(end)=0;

indx=(hx==0);
hx(indx)=[];
indy=(hy==0);
hy(indy)=[];
hx(1)=[];hx(end)=[];
hy(1)=[];hy(end)=[];
hzx=myxcorr(hz-mean(hz),hz-mean(hz));
hxx=myxcorr(hx-mean(hx),hx-mean(hx));
hyx=myxcorr(hy-mean(hy),hy-mean(hy));
ax2=axes(ax.Parent);
subplot(1,2,1,ax);
subplot(1,2,2,ax2);
findx=find(~indx);findy=find(~indy);
plot(ax,zr(2:end-2),hz,zr(2:end-2),hzx/max(hzx)*max(hz));
ax.YLim(2)=(myquantile(hz,.99));
ax.YLim(1)=min(myquantile(hz,.01),myquantile(hzx/max(hzx)*max(hz),.01));
plot(ax2,xr(findx(2:end-1)),hx,xr(findx(2:end-1)),hxx/max(hxx)*max(hx),xr(findy(2:end-1)),hy,xr(findy(2:end-1)),hyx/max(hyx)*max(hy));
end

function zs=testfit(teststack,coeff,p,linepar,ax)
if nargin<4
    linepar={};
elseif ~iscell(linepar)
    linepar={linepar};
end
d=round((size(teststack,1)-p.ROIxy)/2);
            range=d+1:d+p.ROIxy;
            
% if p.doublesample
%     coefffak=4;
%     fitterGPU=@GPUmleFit_LM_v2;
%     fitterCPU=@kernel_MLEfit_Spline_LM_SMAP_v2;
% else
    coefffak=1;
    fitterGPU=@callYimingFitter;
%     fitterCPU=@kernel_MLEfit_Spline_LM_SMAP_v2_nointerp;
% end  
numstack=size(teststack,4);
fprintf('fitting test stacks: 0.00')
    for k=1:size(teststack,4)
        fprintf([ '\b\b\b\b' num2str(k/numstack,'%1.2f')])
%         try
          [P] =  fitterGPU(single(squeeze(teststack(range,range,:,k))),single(coefffak*coeff),100,5,0);
%         catch err
% %             err
%             disp('run on CPU')        
%             [P] =  fitterCPU(teststack(range,range,:,k),(coefffak*coeff),single(p.roisize),100);
%         end
    z=(1:size(P,1))'-1;
    plot(ax,z,P(:,5),linepar{:})
    hold(ax,'on')
    zs(:,k)=P(:,5);
    end
    
end