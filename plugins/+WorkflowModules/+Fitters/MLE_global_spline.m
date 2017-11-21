classdef MLE_global_spline<interfaces.WorkflowFitter
    properties
        fitpar
%         mirrorstack
    end
    methods
        function obj=MLE_global_spline(varargin)
            obj@interfaces.WorkflowFitter(varargin{:})
            obj.inputChannels=1; 
             obj.setInputChannels(1,'frame');
             
        end
        function pard=guidef(obj)
            pard=guidef(obj);
        end
        function fitinit(obj)
            obj.infofields={'x','y','ID','dx','dy'};
            obj.fitpar=getfitpar(obj);
            obj.fitpar.fitfunction=@mleFit_LM_global; %later: include single channel, decide here
            
             transform=obj.getPar('loc_globaltransform');
             
             obj.fitpar.mirrorud=contains(transform.tinfo.mirror.targetmirror,'up');
             obj.fitpar.mirrorrl=contains(transform.tinfo.mirror.targetmirror,'right');
            % check if x,y, then initialize range etc
%             obj.fitpar.fitfunction = @obj.nofound;
%             disp('checking cuda fit')
%             reporttext='GPU fit function did not run. Possibly the wrong CUDA version is installed.';
%             img=zeros(11,'single');img(5,5)=1;
            
%             try
% %                 fitp=callYimingFitter(img,single(1),single(10),single(2),single(0),0);
%                 fitp=mleFit_LM(img,1);
%                 obj.fitpar.fitfunction=@mleFit_LM;
% %                 obj.fitpar.fitfunction=@callYimingFitter;
%                 reporttext='mleFit_LM works';
%             end
            roisize=obj.getPar('loc_ROIsize');
            obj.numberInBlock=round(obj.fitpar.roisperfit*100/roisize^2);
            
            if obj.fitpar.fitmode==5 ||obj.fitpar.fitmode==6
                EMfile=obj.getPar('loc_fileinfo').EMon;
                EMcal=obj.fitpar.EMon;
%                 EMcal=obj.fitpar.splinefit{1}.isEM;
                mirrorstack=obj.getSingleGuiParameter('automirror');
%                 switch mirrorstack.selection
%                     case 'auto'
%                         obj.fitpar.mirrorstack=~(EMfile==EMcal);
%                     case 'mirror'
%                         obj.fitpar.mirrorstack=true;
%                     otherwise
%                         obj.fitpar.mirrorstack=false;
%                 end
%                 if obj.getSingleGuiParameter('automirror')
%                     
%                 else
%                     obj.fitpar.mirrorstack=false;
%                 end
 obj.fitpar.mirrorstack=false; %later: remove completely
                p=obj.getAllParameters;
                if p.overwritePixelsize
                    obj.setPar('overwrite_pixelsize',[p.pixelsizex p.pixelsizey])
                    cs=obj.getPar('loc_cameraSettings');
                    cs.cam_pixelsize_um=[p.pixelsizex p.pixelsizey];
                    obj.setPar('loc_cameraSettings',cs);
                else
                    obj.setPar('overwrite_pixelsize',[])
                end
                 obj.setPar('loc_iterations',p.iterations);
            end   
            
%             disp(reporttext)
        end
        function nofound(obj,varargin)
            disp('fit function not working. Wrong Cuda version?')
        end

        function locs=fit(obj,imstack,stackinfo)
            if obj.fitpar.fitmode==3
                X=stackinfo.X;Y=stackinfo.Y;
                obj.fitpar.zparhere=[obj.fitpar.zpar{X,Y}(:)];
            elseif obj.fitpar.fitmode==5 || obj.fitpar.fitmode==6
                X=stackinfo.X;Y=stackinfo.Y;
                obj.fitpar.splinefithere=[obj.fitpar.splinefit{X,Y}(:)];
            end
            if obj.fitpar.issCMOS
                varstack=getvarmap(obj.fitpar.varmap,stackinfo,size(imstack,1));
            else
                varstack=0;
            end
            out=fitwrapper_global(imstack,obj.fitpar,stackinfo,varstack);
            locs=fit2locs_global(out,stackinfo,obj.fitpar,imstack);
        end

        function initGui(obj)
            initGui@interfaces.WorkflowFitter(obj);
%             obj.guihandles.fitmode.Callback={@fitmode_callback,obj};
%             fitmode_callback(0,0,obj)
            obj.guihandles.loadcal.Callback={@loadcall_callback,obj};
%             obj.addSynchronization('loc_fileinfo',[],[],{@loc_fileinfo_callback,obj});

        end
            
    end
end

function loadscmos_callback(a,b,obj)
fs=obj.getSingleGuiParameter('scmosfile');
if isempty(fs)
    fs='*.*';
end
[file,pfad]=uigetfile(fs);
if file
    obj.setGuiParameters(struct('scmosfile',[pfad file]))
end

end

% function loc_fileinfo_callback(obj)
% 
% end
function locs=fit2locs_global(results,stackinfo,fitpar,image)
if isempty(results)
    locs=[];
    return
end
numl=size(results.P,1);
numpar=5;
numchannels=2;

v1=ones(numl,1,'single');
s=size(image);          
dn=ceil((s(1)-1)/2)*v1;

shiftx=0;%-0.5; %deviation from ground truth
shifty=0;%-0.5;
posx=stackinfo.x(results.indused)+shiftx;
posy=stackinfo.y(results.indused)+shifty;
frame=stackinfo.frame(results.indused);
P=results.P;
EMexcess=fitpar.EMexcessNoise;
CRLB=results.CRLB;
LogL=results.LogL;
           CRLB(isnan(CRLB))= 0; %XXXXXXXXX
           LogL(isnan(LogL))= 0; %XXXXXXXXX
           CRLB((CRLB)<0)= 0; %XXXXXXXXX
          % LogL((LogL)<0)= 0; %XXXXXXXXX
           
o=ones( numpar,1);fac=num2cell(o);z=zeros( numpar,1);off=num2cell(z);
faccrlb=fac;
% locs.xpix=P(:,2)-dn+posx;
if (fitpar.fitmode==5||fitpar.fitmode==6) && fitpar.mirrorstack
    fac{2}=-1;
    off{2}=dn+posx;
%     locs.xpix=dn-P(:,2)+posx;
else
    fac{2}=1;
    off{2}=-dn+posx;
%     locs.xpix=P(:,2)-dn+posx;
end
% locs.ypix=P(:,1)-dn+posy;
off{1}=-dn+posy;


% locs.phot=P(:,4)*EMexcess;
% locs.bg=P(:,5)*EMexcess;
fac{4}=EMexcess;
fac{5}=EMexcess;
faccrlb{4}=EMexcess;
faccrlb{5}=EMexcess;
locs.frame=frame;

% locs.xerrpix=sqrt(CRLB(:,2));
% locs.yerrpix=sqrt(CRLB(:,1));
% locs.photerr=sqrt(CRLB(:,4))*EMexcess;
% locs.bgerr=sqrt(CRLB(:,5))*EMexcess;
locs.logLikelihood=LogL;

locs.peakfindx=posx;
locs.peakfindy=posy;

% switch fitpar.fitmode
%     case 1 %sx not fitted
%         sx=fitpar.PSFx0*v1;
%         locs.PSFxpix=0*locs.xpix+sx;
%         locs.PSFypix=locs.PSFxpix;
%     case 2 % sx free
%         locs.PSFxpix=P(:,5);
%         locs.PSFxerr=sqrt(CRLB(:,5));
% %                     sx=locs.PSFx;
%         locs.PSFypix=locs.PSFxpix;
%     case 3
%         locs.znm=(P(:,5)*1000+fitpar.objPos*v1)*fitpar.refractive_index_mismatch;
%         locs.zerr=sqrt(CRLB(:,5))*1000*fitpar.refractive_index_mismatch;
%         [locs.PSFxpix,locs.PSFypix]=zpar2sigma(locs.znm/1000,fitpar.zparhere);
% 
% 
%     case 4  %sx,sy
% 
%         locs.PSFxpix=P(:,5);
%         locs.PSFxerr=sqrt(CRLB(:,5));
%         locs.PSFypix=P(:,6);
%         locs.PSFyerr=sqrt(CRLB(:,6));  
%     case {5,6}
%         
%         %         locs.znm=(P(:,5)*1000+fitpar.objPos*v1)*fitpar.refractive_index_mismatch;
%         locs.znm=((P(:,3)-fitpar.z0)*fitpar.dz)*fitpar.refractive_index_mismatch;
%         notconverged=P(:,3)<2|P(:,3)>size(fitpar.splinefit{1}.cspline.coeff,3)-2;
%         locs.znm(notconverged)=NaN;
%         
%         locs.zerr=sqrt(CRLB(:,3))*fitpar.dz*fitpar.refractive_index_mismatch;
% %         [locs.PSFxpix,locs.PSFypix]=zpar2sigma(locs.znm/1000,fitpar.zparhere);
%         
         sx=1*v1;
        locs.PSFxpix=sx;
        locs.PSFypix=sx;
% end
fac{3}=fitpar.dz*fitpar.refractive_index_mismatch;
off{3}=-fitpar.z0*fitpar.dz*fitpar.refractive_index_mismatch;
faccrlb{3}=fitpar.dz*fitpar.refractive_index_mismatch;

names={'ypix','xpix','znm','phot','bg'};
linked=fitpar.link;
ind=1;
for k=1:length(names)
    if linked(k)
        locs.(names{k})=P(:,ind).*fac{k}+off{k};
        locs.([names{k} 'err'])=sqrt(CRLB(:,ind)).*faccrlb{k};
        ind=ind+1;
    else
        for c=1:numchannels
            if c==1
                ch='';
            else
                ch=num2str(c);
            end
            locs.([names{k} ch])=P(:,ind).*fac{k}+off{k};
            locs.([names{k} ch 'err'])=sqrt(CRLB(:,ind)).*faccrlb{k};
            ind=ind+1;
        end
    end
end
locs.iterations=P(:,end);
locs.locpthompson=sqrt((locs.PSFxpix.*locs.PSFypix+1/12*v1)./( locs.phot/EMexcess)+8*pi*(locs.PSFxpix.*locs.PSFypix).^2.* locs.bg./( locs.phot/EMexcess).^2);
end

function out=fitwrapper_global(imstack,fitpar,stackinfo,varstack)
numberOfChannels=2;
nfits=size(imstack,3);
npar=5;
numframes=size(imstack,3); 

s=size(imstack);
if length(s)==2 
 s(3)=1;
end
if s(3)<numberOfChannels  %sorting: needs at least two 
    out=[];
 return
end
% fitpar=obj.fitpar;
EMexcess=fitpar.EMexcessNoise;
if isempty(EMexcess)
    EMexcess=1;
end


        
dT=zeros(npar,2,nfits);
dT(1,2,:)=stackinfo.dy;
dT(2,2,:)=stackinfo.dx;
sharedA = repmat(int32(fitpar.link'),[1 numframes]);
 
%for testing: change order


out.indused=1:numberOfChannels:numframes;   

% [~,ind]=sort(rand(length(out.indused),1));


imfit(:,:,:,1)=imstack(:,:,1:numberOfChannels:end);
if fitpar.mirrorud
    imfit(:,:,:,2)=imstack(end:-1:1,:,2:numberOfChannels:end);
    dT(1,2,:)=-dT(1,2,:);
else
    imfit(:,:,:,2)=imstack(:,:,2:numberOfChannels:end);
end

if fitpar.mirrorstack %em vs non em mirror
%     dT(2,2,:)=-dT(2,2,:);
%bit problem with using transform. everything is mirrord.
end

% imfit=imfit(:,:,ind,:);
% dT=dT(:,:,ind);
% dT=zeros(npar,2,nfits);
arguments{1}=imfit;
arguments{2}=sharedA;
arguments{3}=fitpar.iterations;
arguments{4}=single(fitpar.splinefithere.coeff);
arguments{5}=single(dT);
%imstack, sharedflag, iterations, spline coefficients, channelshift,
%fitmode, varmap
arguments{6}=fitpar.fitmode;
[P CRLB LogL]=fitpar.fitfunction(arguments{:});
% htot=P(:,8)
% 
% arguments{5}=varstack;
% arguments{6}=1;
% 
%     switch fitpar.fitmode
%         case {1,2,4} %fix
%             arguments{4}=fitpar.PSFx0;
%             arguments{1}=single(imstack/EMexcess);
% %         case 2 %free
%         case 3 %z
%             arguments{1}=single(imstack/EMexcess);
%             arguments{4}=single(fitpar.zparhere);
% %         case 4 %sx sy
%         case {5,6} %spline   
%             if fitpar.mirrorstack
%                 arguments{1}=single(imstack(:,end:-1:1,:)/EMexcess);
%             else
%                 arguments{1}=single(imstack/EMexcess);
%             end
%             arguments{4}=single(fitpar.splinefithere.cspline.coeff);
%     end
   
%     if fitpar.fitmode==6
%           
%         [P1 CRLB1 LL1 P2 CRLB2 LL2 ]=fitpar.fitfunction(arguments{:});
%         P = repmat(single(LL1>=LL2),1,6).*P1+repmat(single(LL1<LL2),1,6).*P2;
%         CRLB = repmat(single(LL1>=LL2),1,5).*CRLB1+repmat(single(LL1<LL2),1,5).*CRLB2;
%         LogL = repmat(single(LL1>=LL2),1,1).*LL1+repmat(single(LL1<LL2),1,1).*LL2;
%     else
%         [P CRLB LogL]=fitpar.fitfunction(arguments{:});
%     end

out.P=P;
out.CRLB=CRLB;
 out.LogL=LogL;
end
        

function out=fitwrapper(imstack,fitpar,varstack)
s=size(imstack);
if length(s)==2 
 s(3)=1;
end
if s(3)==0
    out=[];
 return
end
% fitpar=obj.fitpar;
EMexcess=fitpar.EMexcessNoise;
if isempty(EMexcess)
    EMexcess=1;
end

arguments{2}=fitpar.fitmode;
arguments{3}=fitpar.iterations;

arguments{5}=varstack;
arguments{6}=1;

    switch fitpar.fitmode
        case {1,2,4} %fix
            arguments{4}=fitpar.PSFx0;
            arguments{1}=single(imstack/EMexcess);
%         case 2 %free
        case 3 %z
            arguments{1}=single(imstack/EMexcess);
            arguments{4}=single(fitpar.zparhere);
%         case 4 %sx sy
        case {5,6} %spline   
            if fitpar.mirrorstack
                arguments{1}=single(imstack(:,end:-1:1,:)/EMexcess);
            else
                arguments{1}=single(imstack/EMexcess);
            end
            arguments{4}=single(fitpar.splinefithere.cspline.coeff);
    end
   
    if fitpar.fitmode==6
          
        [P1 CRLB1 LL1 P2 CRLB2 LL2 ]=fitpar.fitfunction(arguments{:});
        P = repmat(single(LL1>=LL2),1,6).*P1+repmat(single(LL1<LL2),1,6).*P2;
        CRLB = repmat(single(LL1>=LL2),1,5).*CRLB1+repmat(single(LL1<LL2),1,5).*CRLB2;
        LogL = repmat(single(LL1>=LL2),1,1).*LL1+repmat(single(LL1<LL2),1,1).*LL2;
    else
        [P CRLB LogL]=fitpar.fitfunction(arguments{:});
    end

out.P=P;
out.CRLB=CRLB;
out.LogL=LogL;
end
        
        
function loadcall_callback(a,b,obj)
p=obj.getAllParameters;
[f,p]=uigetfile(p.cal_3Dfile);
if f
    l=load([p f]);
    if ~isfield(l,'outforfit') && ~isfield(l,'SXY') && ~isfield(l,'cspline')
        msgbox('no 3D data recognized. Select other file.');
    end
    obj.setGuiParameters(struct('cal_3Dfile',[p f]));
    
end
end

function fitpar=getfitpar(obj)
p=obj.getAllParameters;
fitpar.iterations=p.iterations;
fitpar.fitmode=p.fitmode.Value;
fitpar.roisperfit=p.roisperfit;
fitpar.issCMOS=false;
fitpar.link=p.link;
if fitpar.fitmode==3||fitpar.fitmode==5
    fitpar.issCMOS=p.isscmos;
    fitpar.PSF2D=p.fit2D;
    if p.fit2D
        fitpar.fitmode=6;
    end
    
    calfile=p.cal_3Dfile;
    cal=load(calfile);

    fitpar.objPos=0;
    if isfield(cal,'outforfit')
        fitpar.zpar{1,1}=cal.outforfit;
    elseif isfield(cal,'SXY')
        s=size(cal.SXY);
        Z=1;
%         if p.useObjPos
%             zr=cal.SXY(1).Zrangeall;
%             zr(1)=[];zr(end)=inf;
%             Z=find(p.objPos<=zr,1,'first');
%         end
        for X=s(1):-1:1
            for Y=s(2):-1:1
%                     if isfield(cal.SXY(X,Y,Z),'gauss_zfit')
                zpar{X,Y}=cal.SXY(X,Y,Z).gauss_zfit;
%                     else
%                         zpar{X,Y}=[];
%                     end
                %global:combine splines
%                 cs=cal.SXY(X,Y,Z).cspline_all;
                cs=cal.SXY(X,Y,Z).cspline;
                if iscell(cs.coeff)
                    coeff(:,:,:,:,1)=cs.coeff{1};
                    coeff(:,:,:,:,2)=cs.coeff{2};
                    cs.coeff=single(coeff);
                end
                splinefit{X,Y}=cs;
                

            end
        end
        if ~isempty(splinefit{1})
            fitpar.dz=splinefit{1}.dz;
            fitpar.z0=splinefit{1}.z0;
            fitpar.splinefit=splinefit;
        end
        fitpar.zpar=zpar;

        if numel(cal.SXY)>1
            obj.spatial3Dcal=true;
        else
            obj.spatial3Dcal=false;
        end
        xr=cal.SXY(1,1).Xrangeall;
        xr(1)=-inf;xr(end)=inf;
        yr=cal.SXY(1,1).Yrangeall;
        yr(1)=-inf;yr(end)=inf;
        obj.spatialXrange=xr;
        obj.spatialYrange=yr;
        fitpar.EMon=cal.SXY(1).EMon;
    elseif isfield(cal,'cspline')
        fitpar.zpar{1}=cal.gauss_zfit;
        fitpar.dz=cal.cspline.dz;
        fitpar.z0=cal.cspline.z0;
        fitpar.splinefit{1}=cal.cspline_all;
        if ~isfield(fitpar.splinefit{1}.cspline,'isEM')
            fitpar.splinefit{1}.cspline.isEM=false;
        end
        fitpar.EMon=false;

    else
        disp('no calibration found')

    end

    %load sCMOS
    if p.isscmos
        [~,~,ext]=fileparts(p.scmosfile);
        switch ext
            case '.tif'
                varmaph=imread(p.scmosfile);
            case '.mat'
                varmaph=load(p.scmosfile);
                if isstruct(varmaph)
                    fn=fieldnames(varmaph);
                    varmaph=varmaph.(fn{1});
                end
            otherwise
                disp('could not load variance map. No sCMOS noise model used.')
                p.isscmos=false;
                fitpar.issCMOS=false;
                varstack=0;
                varmaph=[];
        end
        if ~isempty(varmaph)
            roi=p.loc_cameraSettings.roi;
            varmap=varmaph(max(1,roi(1)):roi(3),max(1,roi(2)):roi(4));
        end
    else 
        varmap=[];
    end
    fitpar.varmap=varmap*p.loc_cameraSettings.pix2phot;
    fitpar.refractive_index_mismatch=p.refractive_index_mismatch;


% elseif fitpar.fitmode==5
%     calfile=p.cal_3Dfile;
%     cal=load(calfile);
%     fitpar.splinecoefficients=single(cal.cspline.coeff);
%     fitpar.z0=cal.z0;
%     fitpar.dz=cal.dz; 
%     fitpar.refractive_index_mismatch=p.refractive_index_mismatch;
%     fitpar.objPos=p.objPos;
    
else
    fitpar.PSFx0=p.PSFx0;
end
if p.loc_cameraSettings.EMon
    fitpar.EMexcessNoise=2;
else
fitpar.EMexcessNoise=1;
end

end

function varstack=getvarmap(varmap,stackinfo,roisize)
numim=length(stackinfo.x);
varstack=zeros(roisize,roisize,numim,'single');
dn=floor(roisize/2);
for k=1:numim
%     stackinfo.x(k)
%     stackinfo.y(k)
    varstack(:,:,k)=varmap(stackinfo.x(k)-dn:stackinfo.x(k)+dn,stackinfo.y(k)-dn:stackinfo.y(k)+dn);
end
end

function fitmode_callback(a,b,obj)
p=obj.getGuiParameters;
fitmode=p.fitmode.Value;
% fitz={'loadcal','cal_3Dfile','trefractive_index_mismatch','refractive_index_mismatch','overwritePixelsize','pixelsizex','pixelsizey','automirror','fit2D','isscmos','selectscmos','scmosfile'};
% fitxy={'PSFx0','tPSFx0'};
% switch fitmode
%     case {3,5}
%         ton=fitz;
%         toff=fitxy;
%     otherwise
%         toff=fitz;
%         ton=fitxy;
% end

switch fitmode
    case {1,2}
        roisize=7;
        iterations=30;
      
    otherwise
        roisize=13;
        iterations=50;
end

obj.setPar('loc_ROIsize',roisize);

% obj.fieldvisibility('on',ton,'off',toff);
obj.setGuiParameters(struct('iterations',iterations));
end

function pard=guidef(obj)
p1(1).value=1; p1(1).on={'PSFx0','tPSFx0'}; 
p1(1).off={'loadcal','cal_3Dfile','trefractive_index_mismatch','refractive_index_mismatch','overwritePixelsize',...
    'fit2D','isscmos','pixelsizex','pixelsizey','selectscmos','scmosfile'};
p1(2)=p1(1);p1(2).value=2;
p1(3).value=3;p1(3).off={'PSFx0','tPSFx0'};p1(3).on={'loadcal','cal_3Dfile','trefractive_index_mismatch','refractive_index_mismatch','overwritePixelsize','fit2D','isscmos'};
p1(4)=p1(1);p1(4).value=4;
p1(5)=p1(3);p1(5).value=5;
p1(6)=p1(5);p1(6).value=6;

pard.fitmode.object=struct('Style','popupmenu','String',{{'PSF fix','PSF free','3D z','ellipt: PSFx PSFy','Spline'}},'Value',2,'Callback',{{@obj.switchvisible,p1,{@fitmode_callback,0,0,obj}}});
pard.fitmode.position=[1,1];
pard.fitmode.Width=1.5;
pard.fitmode.TooltipString=sprintf('Fit mode. Fit with constant PSF, free PSF, 3D with astigmatism, asymmetric PSF (for calibrating astigmatic 3D)');

pard.text.object=struct('Style','text','String','Iterations:');
pard.text.position=[1,2.5];
pard.text.Width=0.7;
pard.text.Optional=true;
pard.iterations.object=struct('Style','edit','String','50');
pard.iterations.position=[1,3.2];
pard.iterations.TooltipString=sprintf('number of iterations for the GPU fitter (typical 50, use 100-250 for ellipt: PSFx PSFy or 3Dz).');
pard.iterations.Optional=true;
pard.iterations.Width=0.5;

pard.roisperfitt.object=struct('Style','text','String','ROIs/fit:');
pard.roisperfitt.position=[1,3.9];
pard.roisperfitt.Width=0.6;
pard.roisperfitt.Optional=true;
pard.roisperfit.object=struct('Style','edit','String','15000');
pard.roisperfit.position=[1,4.5];
pard.roisperfit.TooltipString=sprintf('Number of 10 x 10 pixel ROIs passed to GPU for fitting. For other ROI sizes, the number is adjusted accordingly.');
pard.roisperfit.Optional=true;
pard.roisperfitt.TooltipString=pard.roisperfit.TooltipString;
pard.roisperfit.Width=0.5;

pard.tPSFx0.object=struct('Style','text','String','PSFx start (pix)');
pard.tPSFx0.position=[2,1];
pard.tPSFx0.Width=1.25;
pard.tPSFx0.Optional=true;

pard.PSFx0.object=struct('Style','edit','String','1');
pard.PSFx0.position=[2,2.25];
pard.PSFx0.Width=0.75;
pard.PSFx0.TooltipString=sprintf('start value for PSF, or size of PSF when PSF fixed (in camera pixels)');
pard.PSFx0.Optional=true;

pard.fit2D.object=struct('Style','checkbox','String','2D PSF','Value',0);
pard.fit2D.position=[2,3.5];
pard.fit2D.Width=.75;
pard.fit2D.TooltipString=sprintf('Check if PSF model is 2D (no specific PSF engineering), or displays a high degree of similarity above and below the focal plane');
pard.fit2D.Optional=true;

% pard.automirror.object=struct('Style','popupmenu','String',{{'auto','no mirror','mirror'}});
% pard.automirror.position=[2,4.25];
% pard.automirror.Width=.75;
% pard.automirror.Optional=true;

pard.loadcal.object=struct('Style','pushbutton','String','Load 3D cal');
pard.loadcal.position=[2,1];
pard.loadcal.Width=.5;
pard.cal_3Dfile.object=struct('Style','edit','String','settings/cal_3DAcal.mat');
pard.cal_3Dfile.position=[2,1.5];
pard.cal_3Dfile.Width=2;
pard.cal_3Dfile.TooltipString=sprintf('3D calibration file for astigmtic 3D. \n Generate from bead stacks with plugin: Analyze/sr3D/CalibrateAstig');


p(1).value=0; p(1).on={}; p(1).off={'linkt','link'};
p(2).value=1; p(2).on={'linkt','link'}; p(2).off={};

pard.isglobal.object=struct('Style','checkbox','String','Global fit','Callback',{{@obj.switchvisible,p}});
pard.isglobal.position=[3,1];
pard.isglobal.Width=1;
pard.isglobal.Optional=true;

pard.linkt.object=struct('Style','text','String','link: x y z N bg');
pard.linkt.position=[3,2];
pard.linkt.Width=1.5;
pard.linkt.Optional=true;

pard.link.object=struct('Style','edit','String','1 1 1 0 0');
pard.link.position=[3,3.5];
pard.link.Width=1.5;
pard.link.Optional=true;

pard.trefractive_index_mismatch.object=struct('Style','text','String','RI mismatch factor:');
pard.trefractive_index_mismatch.position=[4,3.5];
pard.trefractive_index_mismatch.Width=1.5;
pard.trefractive_index_mismatch.Optional=true;

pard.refractive_index_mismatch.object=struct('Style','edit','String','.8');
pard.refractive_index_mismatch.position=[4,4.5];
pard.refractive_index_mismatch.TooltipString=sprintf('Correction factor to take into account the different refracrive indices of immersion oil and buffer. \n This leads to smaller distances inside the sample compared to bead calibration. \n Bead calibration: in piezo positions (nm). \n This factor transforms z positions to real-space z positions. \n For high-NA oil objectives: typical 0.72 (range 0.7-1).');
pard.refractive_index_mismatch.Optional=true;
pard.refractive_index_mismatch.Width=0.5;


p(1).value=0; p(1).on={}; p(1).off={'pixelsizex','pixelsizey'};
p(2).value=1; p(2).on={'pixelsizex','pixelsizey'}; p(2).off={};
pard.overwritePixelsize.object=struct('Style','checkbox','String','New pixelsize X,Y (um):','Callback',{{@obj.switchvisible,p}});
pard.overwritePixelsize.position=[4,1];
pard.overwritePixelsize.Width=1.5;
pard.overwritePixelsize.Optional=true;

pard.pixelsizex.object=struct('Style','edit','String','.1');
pard.pixelsizex.position=[4,2.5];
pard.pixelsizex.Width=0.5;
pard.pixelsizex.Optional=true;

pard.pixelsizey.object=struct('Style','edit','String','.1');
pard.pixelsizey.position=[4,3];
pard.pixelsizey.Width=0.5;
pard.pixelsizey.Optional=true;

p(1).value=0; p(1).on={}; p(1).off={'selectscmos','scmosfile'};
p(2).value=1; p(2).on={'selectscmos','scmosfile'}; p(2).off={};
pard.isscmos.object=struct('Style','checkbox','String','sCMOS','Callback',{{@obj.switchvisible,p}});   
pard.isscmos.position=[5,1];
pard.isscmos.Optional=true;
pard.selectscmos.object=struct('Style','pushbutton','String','Load var map','Callback',{{@loadscmos_callback,obj}});   
pard.selectscmos.TooltipString='Select sCMOS variance map (in ADU^2) of same size ROI on chip as image stack';
pard.selectscmos.position=[5,2];
pard.selectscmos.Optional=true;
pard.scmosfile.object=struct('Style','edit','String','');
pard.scmosfile.TooltipString='Tiff/.mat image containing sCMOS variance map (same ROI on camera as tiff).';
pard.scmosfile.position=[5,3];
pard.scmosfile.Optional=true;
    pard.scmosfile.Width=2;
    
pard.asymmetry.object=struct('Style','checkbox','String','get asymmetry');   
pard.asymmetry.position=[6,1];
pard.asymmetry.Optional=true;
    
pard.plugininfo.type='WorkflowFitter';
pard.plugininfo.description='Maximum likelyhood estimater, optimized for GPU processing. According to: C. S. Smith, N. Joseph, B. Rieger, and K. A. Lidke, ?Fast, single-molecule localization that achieves theoretically minimum uncertainty.,? Nat Methods, vol. 7, no. 5, pp. 373?375, May 2010.';
end