classdef LocFilter<interfaces.WorkflowModule;
    properties
         pixelsize   
    end
    methods
       function obj=LocFilter(varargin)
            obj@interfaces.WorkflowModule(varargin{:})
            obj.inputChannels=1; 
%              obj.setInputChannels(1,'frame');
        end
        function pard=pardef(obj)
            pard=pardef(obj);
        end
        function initGui(obj)
            initGui@interfaces.WorkflowModule(obj);
        end
        function prerun(obj,p)

           obj.pixelsize=obj.getPar('cameraSettings').pixsize*1000;
            
        end
        function output=run(obj,data,p)
            output=[];
            locs=data.data;%get;
            if isempty(locs)
                disp('no localizations found')
                return
            end
            indin=true(length(locs.frame),1);
            if ~isempty(locs)
                %locprec
                if p.check_locprec && isfield(locs,'xerrpix')
                    phot=locs.xerrpix;
                    if isfield(locs,'yerrpix')
                        phot=sqrt((phot.^2+locs.yerrpix.^2)/2);
                    end
                    xerrnm=phot*obj.pixelsize;
                    val=p.val_locprec;
                    if length(val)==1
                        val=[0 val];
                    end
                    indinh=xerrnm>=val(1)&xerrnm<=val(2);
                    indin=indin&indinh;
                end
                
                %PSFxnm
                if p.check_psf && isfield(locs,'PSFxpix')
                    psf=locs.PSFxpix;
                    if isfield(locs,'PSFypix')
                        psf=sqrt((psf.^2+locs.PSFypix.^2)/2);
                    end
                    psfnm=psf*obj.pixelsize;
                    val=p.val_psf;
                    if length(val)==1
                        val=[0 val];
                    end
                    indinh=psfnm>=val(1)&psfnm<=val(2);
                    indin=indin&indinh;
                end
                
                %phot
                if p.check_phot && isfield(locs,'phot')
                    phot=locs.phot;
                    val=p.val_phot;
                    if length(val)==1
                        val=[val inf];
                    end
                    indinh=phot>=val(1)&phot<=val(2);
                    indin=indin&indinh;
                end 
                
                fn=fieldnames(locs);
                for k=1:length(fn)
                    locsout.(fn{k})=locs.(fn{k})(indin);
                end
                output=data;
                output.data=locsout;
            end
        end

    end
end


function pard=pardef(obj)
pard.check_locprec.object=struct('Style','checkbox','String','xy-locprec (nm)','Value',1);
pard.check_locprec.position=[1,1];
% pard.check1.Width=0.3;
pard.val_locprec.object=struct('Style','edit','String','100');
pard.val_locprec.position=[1,2];
pard.check_psf.object=struct('Style','checkbox','String','PSFxy (nm)','Value',1);
pard.check_psf.position=[2,1];
pard.val_psf.object=struct('Style','edit','String','300');
pard.val_psf.position=[2,2];
pard.check_phot.object=struct('Style','checkbox','String','Photons','Value',0);
pard.check_phot.position=[3,1];
pard.val_phot.object=struct('Style','edit','String','[200 inf]');
pard.val_phot.position=[3,2];
end