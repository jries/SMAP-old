classdef roi2int_expPSF<interfaces.GuiModuleInterface 
    %determines intensity around a localization by a regression of a
    %Gaussian model with fixed positions and sigma. Either amplitude and
    %background or only the amplitude are fitting parameters
    properties
        splinecoeff
        spline
%         transform
        p
        mirror
    end
    methods
        function obj=roi2int_expPSF(varargin)
            obj@interfaces.GuiModuleInterface(varargin{:});
        end
        function pard=guidef(obj)
            pard=guidef(obj);
        end
        function out=evaluate(obj,roi,loc)
            out=roi2int_fit_e(obj,roi,loc);
        end
        function load3Dfile(obj,a,b)
            calf=obj.getSingleGuiParameter('cal_3Dfile');
            if isempty(calf) %nothing selected
                pf=fileparts(obj.getPar('mainfile'));
                calf=[pf,filesep,'*_3dcal.mat'];
            end
            [file,pfad]=uigetfile(calf);
            if file
                obj.setGuiParameters(struct('cal_3Dfile',[pfad file]));
            end
            load_spline(obj)
        end
        function load_spline(obj)
            f=obj.getSingleGuiParameter('cal_3Dfile');
            if isempty(f)
                return
            end
            obj.spline=load(f);
            sppos=obj.guihandles.splinefields.Value;
            sppos=min(length(obj.spline.SXY(sppos).cspline.coeff),sppos); %did not fix calibrator yet..
            obj.splinecoeff=obj.spline.SXY(sppos).cspline.coeff{sppos};
            
        end
        function prerun(obj,p)
             fit3ddir=strrep(pwd,'SMAP','fit3D');
            if exist(fit3ddir,'file')
            addpath(fit3ddir);
            end
            
            fit3ddir=strrep(pwd,'SMAP','fit3Dcspline');
            if exist(fit3ddir,'file')
            addpath(fit3ddir);
            end
            obj.p=obj.getGuiParameters;
            obj.load_spline;
            sppos=obj.guihandles.splinefields.Value;
            transform=obj.spline.transformation;
            if isa(transform,'interfaces.LocTransformN')
                obj.mirror=transformation.mirror(sppos);
            else
                
                obj.mirror=false(1,2);
                 mm=transform.tinfo.mirror.targetmirror;
                 if contains(mm,'up') && sppos==2
                     obj.mirror(2)=true;
                 end
                 if contains(mm,'right') && sppos==2
                      obj.mirror(1)=true;
                 end
            end
%             roi2int_fit_e();
            %TODO include PSF fit
%             global roi2int_fitG_parameters;
%             roi2int_fitG_parameters=obj.getAllParameters;
%             mp=round(sim+1)/2;
%             dn=single(round((roi2int_fitG_parameters.roisize_fit-1)/2));
%             [roi2int_fitG_parameters.X,roi2int_fitG_parameters.Y]=meshgrid(-dn:dn);
            
        end
    end
end


function p=roi2int_fit_e(obj,roi,loc)
% persistent splinehere
% if nargin==0
%     splinehere=[];
% end
% if isempty(splinehere)
% end

%weights not implemented? Do htat!
sim=size(roi);
if length(sim)==2
    sim(3)=1;
end

if obj.mirror(2)
    roi(:,:,:)=roi(end:-1:1,:,:);
    loc.dy=-loc.dy;
end
if obj.mirror(1)
    roi(:,:,:)=roi(:,end:-1:1,:);
    loc.dx=-loc.dx;
end
mp=round(sim+1)/2;
dn=round((obj.p.roisize_fit-1)/2);
p=zeros(sim(3),2,'single');
% n=-dn:dn;
% [X,Y]=meshgrid(n);
% Z=0*X;
% sppos=obj.guihandles.splinefields.Value;
% splinecoeff=obj.spline.SXY(sppos).cspline.coeff{sppos};
zmp=obj.spline.SXY(1).cspline.z0;
    if obj.p.fixz0
        z=zeros(sim(3),1)+obj.p.z0;
    else
        z=loc.z;
    end
    cor=horzcat(loc.dy+dn,loc.dx+dn,z+zmp);
    template = evalSpline(obj.p.roisize_fit,obj.splinecoeff,1,0,cor);
    
for k=1:sim(3)
    templateh=template(:,:,k);
    if obj.p.fitbg
        
        Xmat=horzcat(templateh(:),ones((2*dn+1)^2,1));
        roih=roi(mp(1)-dn:mp(1)+dn,mp(2)-dn:mp(2)+dn,k);
        p(k,:)=Xmat\roih(:);
    else
        Xmat=templateh(:);
        roih=roi(mp(1)-dn:mp(1)+dn,mp(2)-dn:mp(2)+dn,k)-loc.bg(k);
        p(k,1)=Xmat\roih(:);
        p(k,2)=loc.bg(k);
    end
end

end


function pard=guidef(obj)
pard.t1.object=struct('Style','text','String','roisize');
pard.t1.position=[1,1];
pard.t1.TooltipString='Roi size around localizations for fitting';
% pard.t1.Width=0.5;
pard.roisize_fit.object=struct('Style','edit','String','9');
pard.roisize_fit.position=[1,2];
pard.roisize_fit.TooltipString=pard.t1.TooltipString;

pard.t2.object=struct('Style','text','String','spline cal. SXY index:');
pard.t2.position=[2,1];
pard.t2.Width=2;
% pard.t2.TooltipString='Roi size around localizations for fitting';
pard.splinefields.object=struct('Style','popupmenu','String',{{1,2}});
pard.splinefields.position=[2,3];


pard.loadcal_3Dfile.object=struct('Style','pushbutton','String','load','Callback',@obj.load3Dfile);
pard.loadcal_3Dfile.position=[2,4];
pard.loadcal_3Dfile.Width=1;
pard.loadcal_3Dfile.TooltipString=sprintf('3D calibration file for astigmtic 3D. \n Generate from bead stacks with plugin: Analyze/sr3D/CalibrateAstig');


pard.cal_3Dfile.object=struct('Style','edit','String','');
pard.cal_3Dfile.position=[3,1];
pard.cal_3Dfile.Width=4;
pard.cal_3Dfile.TooltipString=sprintf('3D calibration file for astigmtic 3D. \n Generate from bead stacks with plugin: Analyze/sr3D/CalibrateAstig');


pard.fixz0.object=struct('Style','checkbox','String','Fix z pos to frame:');
pard.fixz0.position=[4,1];
pard.fixz0.Width=3;
% pard.fixz0.TooltipString='fix the PSF during the fit (value in nm). Otherwise the fitted PSF size is used.';

pard.z0.object=struct('Style','edit','String','0');
pard.z0.position=[4,4];
% pard.z0.TooltipString=pard.fixpsf.TooltipString;

pard.fitbg.object=struct('Style','checkbox','String','fit BG','Value',1);
pard.fitbg.position=[1,3];
pard.fitbg.Width=2;
% pard.fitbg.TooltipString='If selected, the background is subtracted and the fit is performed with an offset=0. Otherwise the background is a fit parameter.';

info.prefix='fit';
info.name='fit';
info.fields={'fit_n','fit_bg'};
pard.plugininfo=info;
pard.plugininfo.type='WorkflowIntensity'; 
pard.syncParameters={{'cal_3Dfile','cal_3Dfile',{'String'}}};
end
