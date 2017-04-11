classdef ImageFilter<interfaces.WorkflowModule
    properties
        filterkernel
        preview
        filterkernelPSF
    end
    methods
        function obj=ImageFilter(varargin)
            obj@interfaces.WorkflowModule(varargin{:});
            obj.outputParameters={'loc_loc_filter_sigma'};
        end
        function pard=guidef(obj)
            pard=guidef(obj);
        end
%         function initGui(obj)
%             initGui@interfaces.WorkflowModule(obj);
% %            obj.guihandles.loadmetadata.Callback={@loadmetadata_callback,obj};
% %            obj.guihandles.camparbutton.Callback={@camparbutton_callback,obj};
%         end
        function prerun(obj,p)
            p=obj.getAllParameters;
            switch p.filtermode.Value
                case 1 %GAuss
                    fs=p.loc_loc_filter_sigma;
                    if fs>0
                        obj.filterkernel=fspecial('gaussian', max(5,ceil(3.5/2*fs)*2+1), fs);
                    else
                        obj.filterkernel=1;
                    end
                case 2
                    if isempty(obj.filterkernelPSF)
                        error('for PSF filtering you need to load a PSF model _3Dcal.mat before fitting')
                    end
                    obj.filterkernel=obj.filterkernelPSF;
            end
            obj.preview=obj.getPar('loc_preview');
        end
        function dato=run(obj,data,p)
            dato=data;%.copy;
            imf=filter2(obj.filterkernel,data.data);
            if obj.preview&&obj.getPar('loc_previewmode').Value==3&&~isempty(imf)
              
                figure(obj.getPar('loc_outputfig'))
                hold off
                imagesc(imf);
                colorbar;
                axis equal
            end
            dato.data=(imf);
        end
    end
end

function filtermode_callback(object,b,obj)
gauss={'loc_loc_filter_sigma'};
psf={'loadPSF'};
switch object.Value
    case 1 %Gauss
        obj.fieldvisibility('on',gauss,'off',psf);
    case 2 %PSF
        obj.fieldvisibility('off',gauss,'on',psf);
end
end

function loadPSF_callback(object,b,obj)
p=(obj.getPar('lastSMLFile'));
if ~isempty(p)
    p=fileparts(p);
end
[f,p]=uigetfile([p filesep '*.mat']);
if f
    l=load([p f]);
    PSF=l.SXY(1).splinefit.PSF;
    fig=nanmean(PSF,3);
    h=fig-nanmin(fig);
    h=h/nanmax(h(:));
    obj.filterkernelPSF=h;
    
end

end

function pard=guidef(obj)
pard.filtermode.object=struct('Style','popupmenu','String',{{'Gauss: Sigma(pix)','mean PSF'}},'Callback',{{@filtermode_callback,obj}});
pard.filtermode.position=[1,1];
pard.filtermode.Width=1.3;
pard.filtermode.Optional=true;

pard.loc_loc_filter_sigma.object=struct('Style','edit','String','1.2');
pard.loc_loc_filter_sigma.position=[1,2.3];
pard.loc_loc_filter_sigma.Width=.7;
pard.loc_loc_filter_sigma.TooltipString=sprintf('Sigma (in camera pixels) for a Gaussian filter which is applied after background correction and before peak finding. \n Typical size of PSF in pixels, eg 1 (range: 0.5-5) ');

pard.loadPSF.object=struct('Style','pushbutton','String','load','Callback',{{@loadPSF_callback,obj}},'Visible','off');
pard.loadPSF.position=[1,2.3];
pard.loadPSF.Width=.7;

pard.plugininfo.type='WorkflowModule';
pard.loc_loc_filter_sigma.Optional=true;
pard.plugininfo.description='Gaussian filter, usually applied after background correction and before peak finding.';
pard.text.TooltipString=pard.loc_loc_filter_sigma.TooltipString;
end