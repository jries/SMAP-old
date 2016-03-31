classdef CameraConverter<interfaces.WorkflowModule
    properties
        calfile='settings/CameraCalibration.xls';
        cameraSettings=struct('camId','default','port','Conventional','exposure',1,'emgain',1,'conversion',1,'offset',400,'pixsize',0.1,...
            'roi',[],'temperature',0,'timediff',0);
        cameraSettingsStructure=struct('camId','default','port','Conventional','exposure',1,'emgain',1,'conversion',1,'offset',400,'pixsize',0.1,...
            'roi',[],'temperature',0,'timediff',0);
        EMexcessNoise;
    end
    methods
        function obj=CameraConverter(varargin)
            obj@interfaces.WorkflowModule(varargin{:})
            obj.propertiesToSave={'cameraSettings','EMexcessNoise'};
        end
        function pard=pardef(obj)
            pard=pardef;
        end
        function initGui(obj)
           initGui@interfaces.WorkflowModule(obj);
           obj.guihandles.loadmetadata.Callback={@loadmetadata_callback,obj};
           obj.guihandles.camparbutton.Callback={@camparbutton_callback,obj};
           obj.outputParameters={'cameraSettings','EMexcessNoise'};
           obj.addSynchronization('loc_metadatafile',obj.guihandles.metadatafile,'String',@obj.readmetadata)
%            str=obj.guihandles.metadatafile.String;
%            if exist(str,'file')
%                obj.readmetadata;
%            end
        end
        function prerun(obj,p)
           
        end
        function datao=run(obj,data,p)
%             outpu?t=[];
            pc=obj.cameraSettings;
            offset=pc.offset;
            adu2phot=(pc.conversion/pc.emgain);
            imphot=single((data.data-offset))*adu2phot;
%             imphot=single((data.get-offset))*adu2phot;
%             datao=data.copy;
               datao=data;
               datao.data=imphot;
%             figure(88);imagesc(imphot);
%             datao.set(imphot);
%             drawnow
%             obj.output(dat?ao);
            
            
        end
        
        function readmetadata(obj)
             p=obj.getGuiParameters;
            camcalib=readtable(obj.calfile);
            fid=fopen(p.metadatafile);
            if fid>0
                minfo=fread(fid,[1,100000],'*char');
                fclose(fid);
                metadata=minfoparsec(minfo,camcalib);
                obj.status('metadata updated');drawnow;
            else 
                obj.status('no metadata found, using previous metadata. ');drawnow;
                metadata=[];
            end
            if ~isempty(metadata) %otherwise: keep default
                obj.cameraSettings=copyfields(obj.cameraSettings,metadata,fieldnames(obj.cameraSettings));
            else
                %default settings
            end
            
            if strcmp(obj.cameraSettings.port,'Conventional')|| strcmp(obj.cameraSettings.port,'Normal') %still valid?
                emexcess=1;
            else
                emexcess=2;
            end
            obj.EMexcessNoise=emexcess;
            obj.updateObjectParameters;
        end
    end
end



function camparbutton_callback(a,b,obj)
fn=fieldnames(obj.cameraSettingsStructure);
for k=1:length(fn)
    fields{k}=fn{k};
    defAns{k}=num2str(obj.cameraSettings.(fn{k}));
end
answer=inputdlg(fields,'Acquisition settings',1,defAns);
if ~isempty(answer)
    for k=1:length(fn)
        if isnumeric(obj.cameraSettings.(fn{k}))
            obj.cameraSettings.(fn{k})=str2num(answer{k});
        else
            obj.cameraSettings.(fn{k})=(answer{k});
        end
    end
%     obj.globpar.parameters.cameraSettings=obj.cameraSettings; %doesnt work
end
end

function loadmetadata_callback(a,b,obj)
[f,p]=uigetfile('*.*','Select metadata.txt, tif or _sml.mat');
[~,~,ext ]=fileparts(f);
switch ext
    case '.mat'
        
    case '.tif'
        imloader=imageLoader([p f]);
        metafile=imloader.info.metafile;
        par.metadatafile=metafile;
        obj.setGuiParameters(par);
        obj.updateGui;
    case '.txt'
        par.metadatafile=[p f];
        obj.setGuiParameters(par);
        obj.updateGui;
end
end

% function parseMetatdata(obj)
% fid=fopen(file);
% metadatatxt=fread(fid,[1,10000],'*char');
% fclose(fid);
% minfo=minfoparsec(metadatatxt);
% end


function pard=pardef

pard.text.object=struct('Style','text','String','Acquisition parameters');
pard.text.position=[1,1];
pard.text.Width=3;

pard.metadatafile.object=struct('Style','edit','String',' ');
pard.metadatafile.position=[2,1];
pard.metadatafile.Width=4;

pard.loadmetadata.object=struct('Style','pushbutton','String','Load metadata');
pard.loadmetadata.position=[3,1];

pard.camparbutton.object=struct('Style','pushbutton','String','set Cam Parameters');
pard.camparbutton.position=[3,3.5];
pard.camparbutton.Width=1.5;
end