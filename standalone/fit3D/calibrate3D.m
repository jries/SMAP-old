function calibrate3D(p)
% p.filelist
% p.outputfile
% p.dz
% p.modality
% p.zcorr
% p.ROIxy
% p.ROIz
% p.smoothxy
% p.smoothz
% p.gaussrange
% p.filter;
% p.zcorrframes
% p.gaussroi

%get bead positions
p.status.String='load files and segmet beads';drawnow
f=figure('Name','Bead calibration');
p.tabgroup=uitabgroup(f);
%get beads from images
[beads,p]=images2beads_so(p);
p.midpoint=round(size(beads(1).stack.image,3)/2); %reference for beads
p.ploton=false;

% beads=beads([2 3])
% beads(2)=beads(1);
% beads(3:end)=[];

if contains(p.modality,'astig') || contains(p.modality,'2D')
    %determine sx,sy
%     disp('fit beads to get sx,sy')
    t=tic;
    p.status.String=['Gaussian fit of beads to get spatial paramters '];drawnow
    for k=1:length(beads)
        stackh=single(beads(k).stack.image);
        s=size(stackh); 
        d=round((s(1)-p.gaussroi)/2);
        stack=stackh(d+1:end-d,d+1:end-d,:);
        %fit bead bead stacks with Gaussian model
        if contains(p.modality,'astig')
            P=mleFit_LM(stack,4,100,1,0,1);
            beads(k).loc.PSFxpix=P(:,5);
            beads(k).loc.PSFypix=P(:,6);
            beads(k).loc.phot=P(:,3);
            beads(k).f0=stackas2z_so(beads(k).loc.PSFxpix,beads(k).loc.PSFypix,beads(k).loc.frames,beads(k).loc.phot,p);
        else
            P=mleFit_LM(stack,2,100,1,0,1);
            beads(k).loc.PSFxpix=P(:,5);
            beads(k).loc.PSFypix=P(:,5);
            beads(k).loc.phot=P(:,3);
            beads(k).f0=stackas2z2D_so(beads(k).loc.PSFxpix,beads(k).loc.frames,beads(k).loc.phot,p);
        end
        
        beads(k).loc.bg=P(:,4);
        %determine true position of the beads as the position, where PSFx==PSFy
        
        ind=find(beads(k).loc.frames<=beads(k).f0,1,'last');
        if isnan(beads(k).f0)||isempty(ind)
            ind=1;
        end
        beads(k).psfx0=beads(k).loc.PSFxpix(ind);
        beads(k).psfy0=beads(k).loc.PSFypix(ind);
        if toc(t)>1
            p.status.String=['Gaussian fit of beads to get spatial paramters: ' num2str(k) ' of ' num2str(length(beads))];
            drawnow
            t=tic;
        end
    end
    %remove beads for which no position could be found
    badind=isnan([beads(:).f0]);
    beads(badind)=[];
else
    f0g=p.midpoint;
    for k=1:length(beads)
        beads(k).f0=f0g;
    end
end
if contains(p.modality,'astig')
    %get calibration for Gaussian fit
    p.status.String='get spline approximation';drawnow
    p.ax_z=axes(uitab(p.tabgroup,'Title','sx(z), sy(z)'));
    [spline_curves,indgoodc,curves]=getspline_so(beads,p); 
    gausscal.spline_curves=spline_curves;
    drawnow
else
    indgoodc=true(size(beads));
    gausscal=[];
end
% get cspline calibration
p.status.String='get cspline calibration';drawnow
[csplinecal,indgoods]=getstackcal_so(beads(indgoodc),p);
icf=find(indgoodc);
icfs=icf(indgoods);
cspline.coeff=single(csplinecal.cspline.coeff);
cspline.dz=csplinecal.cspline.dz;
cspline.z0=csplinecal.cspline.z0;
cspline.x0=csplinecal.cspline.x0;

if contains(p.modality,'astig')
    photbead=10^5; %corr PSF normalized to 1. As MLE is used, this screws up statistics totally. Thus assign bright signal to bead.
    stackb=csplinecal.PSF;
    stackb=(stackb)*photbead;
    mp=ceil(size(stackb,1)/2);dx=floor(p.gaussroi/2);
    
    stack=single(stackb(mp-dx:mp+dx,mp-dx:mp+dx,:));
    P=mleFit_LM(stack,4,200,1,0,1);
    ch.sx=double(P(:,5));
    ch.sy=double(P(:,6));
    f0m=median([beads(icfs).f0]);
    ch.z=double(((1:size(stack,3))'-f0m)*p.dz);
    
    p.ax_sxsy=axes(uitab(p.tabgroup,'Title','sx^2-sy^2'));
    p.ax_z.NextPlot='add';
    p.status.String='get Gauss model calibration';drawnow
    gausscalh=getgausscal_so(ch,p); 
    legend(p.ax_z,'bad bead data','good bead data','spline fit sx','spline fit sy','average PSF','average PSF','Gauss zfit','Gauss zfit')

    gausscal=copyfields(gausscal,gausscalh);
    gauss_zfit=single(gausscal.fitzpar);
    gauss_sx2_sy2=gausscal.Sx2_Sy2;
else
    gausscal=[];
    gauss_sx2_sy2=[];
    gauss_zfit=[];
end
cspline_all=csplinecal;
p.status.String='save calibration';drawnow
save(p.outputfile,'gausscal','cspline_all','gauss_sx2_sy2','gauss_zfit','cspline');
p.status.String='Calibration done';drawnow
end




