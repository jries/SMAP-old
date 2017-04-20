classdef MLE_GPU_Yiming<interfaces.WorkflowFitter
    properties
        fitpar
%         mirrorstack
    end
    methods
        function obj=MLE_GPU_Yiming(varargin)
            obj@interfaces.WorkflowFitter(varargin{:})
            obj.inputChannels=1; 
             obj.setInputChannels(1,'frame');
        end
        function pard=guidef(obj)
            pard=guidef;
        end
        function fitinit(obj)
            obj.fitpar=getfitpar(obj);
            % check if x,y, then initialize range etc
            obj.fitpar.fitfunction = @obj.nofound;
            disp('checking cuda fit')
            reporttext='GPU fit function did not run. Possibly the wrong CUDA version is installed.';
            img=zeros(11,'single');img(5,5)=1;
            
            try
                fitp=callYimingFitter(img,single(1),single(10),single(2),single(0),0);
                obj.fitpar.fitfunction=@callYimingFitter;
%                 obj.fitpar.fitfunction=@callYimingFitter;
                 reporttext='GPUmleFit_LM works';
            end
            roisize=obj.getPar('loc_ROIsize');
            obj.numberInBlock=round(obj.fitpar.roisperfit*100/roisize^2);
            
            if obj.fitpar.fitmode==5
                EMfile=obj.getPar('loc_fileinfo').EMon;
                EMcal=obj.fitpar.splinefit{1}.cspline.isEM;
                obj.fitpar.mirrorstack=~(EMfile==EMcal);
                p=obj.getAllParameters;
                if p.overwritePixelsize
                    obj.setPar('overwrite_pixelsize',[p.pixelsizex p.pixelsizey])
                    cs=obj.getPar('loc_cameraSettings');
                    cs.cam_pixelsize_um=[p.pixelsizex p.pixelsizey];
                    obj.setPar('loc_cameraSettings',cs);
                else
                    obj.setPar('overwrite_pixelsize',[])
                end
            end   
            
            disp(reporttext)
        end
        function nofound(obj,varargin)
            disp('fit function not working. Wrong Cuda version?')
        end

        function locs=fit(obj,imstack,stackinfo)
            if obj.fitpar.fitmode==3
                X=stackinfo.X;Y=stackinfo.Y;
                obj.fitpar.zparhere=[obj.fitpar.zpar{X,Y}(:)];
            elseif obj.fitpar.fitmode==5
                X=stackinfo.X;Y=stackinfo.Y;
                obj.fitpar.splinefithere=[obj.fitpar.splinefit{X,Y}(:)];
            end
            out=fitwrapper(imstack,obj.fitpar);
            locs=fit2locs(out,stackinfo,obj.fitpar,imstack);
        end

        function initGui(obj)
            initGui@interfaces.WorkflowFitter(obj);
            obj.guihandles.fitmode.Callback={@fitmode_callback,obj};
            fitmode_callback(0,0,obj)
            obj.guihandles.loadcal.Callback={@loadcall_callback,obj};
%             obj.addSynchronization('loc_fileinfo',[],[],{@loc_fileinfo_callback,obj});

        end
            
    end
end

% function loc_fileinfo_callback(obj)
% 
% end
function locs=fit2locs(results,stackinfo,fitpar,image)
if isempty(results)
    locs=[];
    return
end
numl=size(results.P,1);

v1=ones(numl,1,'single');
s=size(image);          
dn=ceil((s(1)-1)/2)*v1;

shiftx=0;%-0.5; %deviation from ground truth
shifty=0;%-0.5;
posx=stackinfo.x+shiftx;
posy=stackinfo.y+shifty;
frame=stackinfo.frame;
P=results.P;
EMexcess=fitpar.EMexcessNoise;
CRLB=results.CRLB;
LogL=results.LogL;
           CRLB(isnan(CRLB))= 0; %XXXXXXXXX
           LogL(isnan(LogL))= 0; %XXXXXXXXX
           CRLB((CRLB)<0)= 0; %XXXXXXXXX
           LogL((LogL)<0)= 0; %XXXXXXXXX
           
           
% locs.xpix=P(:,2)-dn+posx;
if fitpar.fitmode==5&& fitpar.mirrorstack
    locs.xpix=dn-P(:,2)+1+posx;
else
    locs.xpix=P(:,2)-dn+posx;
end
locs.ypix=P(:,1)-dn+posy;

locs.phot=P(:,3)*EMexcess;
locs.bg=P(:,4)*EMexcess;
locs.frame=frame;

locs.xerrpix=sqrt(CRLB(:,2));
locs.yerrpix=sqrt(CRLB(:,1));
locs.photerr=sqrt(CRLB(:,3))*EMexcess;
locs.bgerr=sqrt(CRLB(:,4))*EMexcess;
locs.logLikelihood=LogL;

locs.peakfindx=posx;
locs.peakfindy=posy;

switch fitpar.fitmode
    case 1 %sx not fitted
        sx=fitpar.PSFx0*v1;
        locs.PSFxpix=0*locs.xpix+sx;
        locs.PSFypix=locs.PSFxpix;
    case 2 % sx free
        locs.PSFxpix=P(:,5);
        locs.PSFxerr=sqrt(CRLB(:,5));
%                     sx=locs.PSFx;
        locs.PSFypix=locs.PSFxpix;
    case 3
        locs.znm=(P(:,5)*1000+fitpar.objPos*v1)*fitpar.refractive_index_mismatch;
        locs.zerr=sqrt(CRLB(:,5))*1000*fitpar.refractive_index_mismatch;
        [locs.PSFxpix,locs.PSFypix]=zpar2sigma(locs.znm/1000,fitpar.zparhere);


    case 4  %sx,sy

        locs.PSFxpix=P(:,5);
        locs.PSFxerr=sqrt(CRLB(:,5));
        locs.PSFypix=P(:,6);
        locs.PSFyerr=sqrt(CRLB(:,6));  
    case 5
        
        %         locs.znm=(P(:,5)*1000+fitpar.objPos*v1)*fitpar.refractive_index_mismatch;
        locs.znm=((P(:,5)-fitpar.z0)*fitpar.dz)*fitpar.refractive_index_mismatch;
        notconverged=P(:,5)<2|P(:,5)>size(fitpar.splinefit{1}.cspline.coeff,3)-2;
        locs.znm(notconverged)=NaN;
        
        locs.zerr=sqrt(CRLB(:,5))*fitpar.dz*fitpar.refractive_index_mismatch;
%         [locs.PSFxpix,locs.PSFypix]=zpar2sigma(locs.znm/1000,fitpar.zparhere);
        
         sx=1*v1;
        locs.PSFxpix=sx;
        locs.PSFypix=sx;
end
locs.locpthompson=sqrt((locs.PSFxpix.*locs.PSFypix+1/12*v1)./( locs.phot/EMexcess)+8*pi*(locs.PSFxpix.*locs.PSFypix).^2.* locs.bg./( locs.phot/EMexcess).^2);
end

function out=fitwrapper(imstack,fitpar)


%             if any(imstack(:)<0)
%                 figure(88);
%                 imageslicer(imstack);
% 
%                 
%             end
            
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

arguments{3}=fitpar.iterations;
arguments{4}=fitpar.fitmode;
arguments{5}=fitpar.issCMOS;
arguments{6}=1;
% try   
    switch fitpar.fitmode
        case {1,2,4} %fix
            arguments{2}=fitpar.PSFx0;
            arguments{1}=single(imstack/EMexcess);
%         case 2 %free
        case 3 %z
            arguments{1}=single(imstack/EMexcess);
            arguments{2}=fitpar.zparhere(1);
            arguments(7:13)=num2cell(fitpar.zparhere(2:8));
%         case 4 %sx sy
        case 5 %spline   
            if fitpar.mirrorstack
                arguments{1}=single(imstack(:,end:-1:1,:)/EMexcess);
            else
                arguments{1}=single(imstack/EMexcess);
            end
            arguments{2}=single(fitpar.splinefithere.cspline.coeff);
    end
    
    [P CRLB LogL]=fitpar.fitfunction(arguments{:});
   

out.P=P;
out.CRLB=CRLB;
out.LogL=LogL;
end
        
        
function loadcall_callback(a,b,obj)
p=obj.getAllParameters;
[f,p]=uigetfile(p.cal_3Dfile);
if f
    l=load([p f]);
    if ~isfield(l,'outforfit') && ~isfield(l,'SXY')
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
if fitpar.fitmode==3||fitpar.fitmode==5
    calfile=p.cal_3Dfile;
    cal=load(calfile);
    if 0% p.useObjPos
        
        disp('obj. position not implemented yet')
    else
        fitpar.objPos=p.objPos;
        if isfield(cal,'outforfit')
            fitpar.zpar{1,1}=cal.outforfit;
        else
            s=size(cal.SXY);
            Z=1;
            if p.useObjPos
                zr=cal.SXY(1).Zrangeall;
                zr(1)=[];zr(end)=inf;
                Z=find(p.objPos<=zr,1,'first');
            end
            for X=s(1):-1:1
                for Y=s(2):-1:1
                    zpar{X,Y}=cal.SXY(X,Y,Z).fitzpar;
                    splinefit{X,Y}=cal.SXY(X,Y,Z).splinefit;
                    
                end
            end
            if ~isempty(splinefit{1})
                fitpar.dz=splinefit{1}.cspline.dz;
                fitpar.z0=splinefit{1}.cspline.z0;
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
                
        end
        fitpar.refractive_index_mismatch=p.refractive_index_mismatch;
    end
    
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
fitpar.issCMOS=false;
end

function fitmode_callback(a,b,obj)
p=obj.getGuiParameters;
fitmode=p.fitmode.Value;
fitz={'loadcal','cal_3Dfile','useObjPos','objPos','trefractive_index_mismatch','refractive_index_mismatch','overwritePixelsize','pixelsizex','pixelsizey'};
fitxy={'PSFx0','tPSFx0'};
switch fitmode
    case {3,5}
        ton=fitz;
        toff=fitxy;
    otherwise
        toff=fitz;
        ton=fitxy;
end

switch fitmode
    case {1,2}
        roisize=7;
        iterations=50;
      
    otherwise
        roisize=13;
        iterations=50;
end

obj.setPar('loc_ROIsize',roisize);

obj.fieldvisibility('on',ton,'off',toff);
obj.setGuiParameters(struct('iterations',iterations));
end

function pard=guidef
pard.fitmode.object=struct('Style','popupmenu','String',{{'PSF fix','PSF free','3D z','ellipt: PSFx PSFy','Spline'}},'Value',2);
pard.fitmode.position=[1,1];
pard.fitmode.Width=2;
pard.fitmode.TooltipString=sprintf('Fit mode. Fit with constant PSF, free PSF, 3D with astigmatism, asymmetric PSF (for calibrating astigmatic 3D)');

pard.text.object=struct('Style','text','String','Iterations:');
pard.text.position=[1,3.3];
pard.text.Width=0.7;
pard.text.Optional=true;
pard.iterations.object=struct('Style','edit','String','50');
pard.iterations.position=[1,4];
pard.iterations.TooltipString=sprintf('number of iterations for the GPU fitter (typical 50, use 100-250 for ellipt: PSFx PSFy or 3Dz).');
pard.iterations.Optional=true;

pard.roisperfitt.object=struct('Style','text','String','ROIs per fit:');
pard.roisperfitt.position=[2,3.3];
pard.roisperfitt.Width=0.7;
pard.roisperfitt.Optional=true;
pard.roisperfit.object=struct('Style','edit','String','5000');
pard.roisperfit.position=[2,4];
pard.roisperfit.TooltipString=sprintf('Number of 10 x 10 pixel ROIs passed to GPU for fitting. For other ROI sizes, the number is adjusted accordingly.');
pard.roisperfit.Optional=true;
pard.roisperfitt.TooltipString=pard.roisperfit.TooltipString;

pard.tPSFx0.object=struct('Style','text','String','PSFx start (pix)');
pard.tPSFx0.position=[2,1];
pard.tPSFx0.Width=1.25;
pard.tPSFx0.Optional=true;

pard.PSFx0.object=struct('Style','edit','String','1');
pard.PSFx0.position=[2,2.25];
pard.PSFx0.Width=0.75;
pard.PSFx0.TooltipString=sprintf('start value for PSF, or size of PSF when PSF fixed (in camera pixels)');
pard.PSFx0.Optional=true;

pard.loadcal.object=struct('Style','pushbutton','String','Load 3D cal');
pard.loadcal.position=[3,1];
pard.cal_3Dfile.object=struct('Style','edit','String','settings/cal_3DAcal.mat');
pard.cal_3Dfile.position=[3,2];
pard.cal_3Dfile.Width=3;
pard.cal_3Dfile.TooltipString=sprintf('3D calibration file for astigmtic 3D. \n Generate from bead stacks with plugin: Analyze/sr3D/CalibrateAstig');

pard.useObjPos.object=struct('Style','checkbox','String','Use objective position:');
pard.useObjPos.position=[4,1];
pard.useObjPos.Width=1.5;
pard.useObjPos.Optional=true;

pard.objPos.object=struct('Style','edit','String','0');
pard.objPos.position=[4,2.5];
pard.objPos.TooltipString=sprintf('Position of the objective above the coverslip (nm, piezo position). \n Only used in combination with CalibrateAstigDeep.');
pard.objPos.Optional=true;
pard.objPos.Width=0.5;

pard.trefractive_index_mismatch.object=struct('Style','text','String','RI mismatch factor:');
pard.trefractive_index_mismatch.position=[4,3];
pard.trefractive_index_mismatch.Width=1.5;
pard.trefractive_index_mismatch.Optional=true;

pard.refractive_index_mismatch.object=struct('Style','edit','String','.8');
pard.refractive_index_mismatch.position=[4,4.5];
pard.refractive_index_mismatch.TooltipString=sprintf('Correction factor to take into account the different refracrive indices of immersion oil and buffer. \n This leads to smaller distances inside the sample compared to bead calibration. \n Bead calibration: in piezo positions (nm). \n This factor transforms z positions to real-space z positions. \n For high-NA oil objectives: typical 0.72 (range 0.7-1).');
pard.refractive_index_mismatch.Optional=true;
pard.refractive_index_mismatch.Width=0.5;

pard.overwritePixelsize.object=struct('Style','checkbox','String','Overwrite pixelsize X,Y (um):');
pard.overwritePixelsize.position=[5,1];
pard.overwritePixelsize.Width=1.5;
pard.overwritePixelsize.Optional=true;

pard.pixelsizex.object=struct('Style','edit','String','.1');
pard.pixelsizex.position=[5,2.5];
pard.pixelsizex.Width=0.5;
pard.pixelsizex.Optional=true;

pard.pixelsizey.object=struct('Style','edit','String','.1');
pard.pixelsizey.position=[5,3];
pard.pixelsizey.Width=0.5;
pard.pixelsizey.Optional=true;


pard.plugininfo.type='WorkflowFitter';
pard.plugininfo.description='Maximum likelyhood estimater, optimized for GPU processing. According to: C. S. Smith, N. Joseph, B. Rieger, and K. A. Lidke, ?Fast, single-molecule localization that achieves theoretically minimum uncertainty.,? Nat Methods, vol. 7, no. 5, pp. 373?375, May 2010.';
end