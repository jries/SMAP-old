classdef CameraConverter<interfaces.WorkflowModule
    properties
        calfile='settings/CameraCalibration.xls';
        loc_cameraSettings=interfaces.metadataSMAP;
%         struct('camId','default','port','Conventional','exposure',1,'emgain',1,'conversion',1,'offset',400,'pixsize',0.1,...
%             'roi',[],'temperature',0,'timediff',0,'comment','');
        loc_cameraSettingsStructure=struct('EMon',1,'emgain',1,'conversion',1,'offset',400,'pixsize',0.1,...
            'roi',[],'exposure',1,'timediff',0,'comment','');
        EMexcessNoise;
        calibrategain=false;
        calibratecounter
        calibrateimages
    end
    methods
        function obj=CameraConverter(varargin)
            obj@interfaces.WorkflowModule(varargin{:})
            obj.loc_cameraSettings=interfaces.metadataSMAP;
            obj.propertiesToSave={'loc_cameraSettings'};
        end
        function pard=guidef(obj)
            pard=guidef;
        end
        function initGui(obj)
           initGui@interfaces.WorkflowModule(obj);
           obj.guihandles.loadmetadata.Callback={@loadmetadata_callback,obj};
           obj.guihandles.camparbutton.Callback={@camparbutton_callback,obj};
           obj.guihandles.calibrate.Callback={@calibrate_callback,obj};
            obj.outputParameters={'loc_cameraSettings'};
           obj.addSynchronization('loc_fileinfo',[],[],@obj.setmetadata)
        end
        function setmetadata(obj)
            md=obj.loc_cameraSettings;
            obj.loc_cameraSettings=interfaces.metadataSMAP;
            obj.loc_cameraSettings=copyfields(obj.loc_cameraSettings,md);
            settings=obj.getPar('loc_fileinfo');
            fn=fieldnames(settings);
            for k=1:length(fn)
                if settings.assigned.(fn{k})
                    obj.loc_cameraSettings.(fn{k})=settings.(fn{k});
                end
            end
            obj.setPar('loc_cameraSettings',obj.loc_cameraSettings);
            obj.setPar('EMon',obj.loc_cameraSettings.EMon);
%             if obj.loc_cameraSettings.EMon
%             	obj.EMexcessNoise=2;
%             else
%                 obj.EMexcessNoise=1;
%             end
        end
        function prerun(obj,p)
           
        end
        function datao=run(obj,data,p)
            
            %calibrate gain offset from images
            global SMAP_stopnow           
             
            if obj.calibrategain
                datao=[];
                if data.eof            
                    SMAP_stopnow=false;
                    calculategain(obj.calibrateimages);
                    obj.calibrategain=false;
                    return
                end

                obj.calibrateimages(:,:,obj.calibratecounter)=data.data;
                obj.calibratecounter=obj.calibratecounter+1;
                if obj.calibratecounter>140
                    SMAP_stopnow=true;
                end
                
                
                
                return
            end
            pc=obj.loc_cameraSettings;
            offset=pc.offset;
            adu2phot=(pc.conversion/pc.emgain);
            imphot=single((data.data-offset))*adu2phot;
               datao=data;
               datao.data=imphot;  
        end
        
        function readmetadata(obj,file)
            dagl
             p=obj.getGuiParameters;
             p.metadatafile=fiele;
%             if ~isfield(p,'metadatafile')
%                 return
%             end
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
                obj.loc_cameraSettings=copyfields(obj.loc_cameraSettings,metadata);%,fieldnames(obj.loc_cameraSettings));
            else
                %default settings
            end
            
            if strcmp(obj.loc_cameraSettings.port,'Conventional')|| strcmp(obj.loc_cameraSettings.port,'Normal') %still valid?
                emexcess=1;
            else
                emexcess=2;
            end
            obj.EMexcessNoise=emexcess;
            obj.updateObjectParameters;
        end
    end
end

function calibrate_callback(a,b,obj)
obj.calibrategain=true;
obj.calibratecounter=1;
obj.calibrateimages=[];
obj.parent.run;
end

function camparbutton_callback(a,b,obj)
fn=fieldnames(obj.loc_cameraSettingsStructure);
%remove later: only there because parameters saved with workflow don't
%include comments
% if ~isfield(obj.loc_cameraSettings,'comment')
%     obj.loc_cameraSettings.comment='no comments';
% end
%XXX
for k=1:length(fn)
    fields{k}=fn{k};
    defAns{k}=num2str(obj.loc_cameraSettings.(fn{k}));
end
answer=inputdlg(fields,'Acquisition settings',1,defAns);
if ~isempty(answer)
    for k=1:length(fn)
        if isnumeric(obj.loc_cameraSettings.(fn{k}))||islogical(obj.loc_cameraSettings.(fn{k}))
            obj.loc_cameraSettings.(fn{k})=str2num(answer{k});
        else
            obj.loc_cameraSettings.(fn{k})=(answer{k});
        end
    end
%     obj.setPar('loc_cameraSettings',obj.loc_cameraSettings);
obj.setmetadata;
%     if obj.loc_cameraSettings.EMon
%         obj.EMexcessNoise=2;
%     else
%         obj.EMexcessNoise=1;
%     end
%     obj.globpar.parameters.loc_cameraSettings=obj.loc_cameraSettings; %doesnt work
end
end

function loadmetadata_callback(a,b,obj)
disp('under construction')
return
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
obj.readmetadata;
end

% function parseMetatdata(obj)
% fid=fopen(file);
% metadatatxt=fread(fid,[1,10000],'*char');
% fclose(fid);
% minfo=minfoparsec(metadatatxt);
% end
function calculategain(img)
img=double(img);
m=mean(img,3);
v=var(img,0,3);
figure(88);
plot(m(:),v(:),'.')
gain=15.6;offs=200;em=200;
pix2phot=gain/em;
hold on
plot(m(:),1/pix2phot*(m(:)-offs),'.')
hold off
end

function pard=guidef

pard.text.object=struct('Style','text','String','Acquisition parameters');
pard.text.position=[1,1];
pard.text.Width=3;
pard.text.Optional=true;
% pard.metadatafile.object=struct('Style','edit','String',' ');
% pard.metadatafile.position=[2,1];
% pard.metadatafile.Width=4;
% pard.metadatafile.Optional=true;
pard.loadmetadata.object=struct('Style','pushbutton','String','Load metadata');
pard.loadmetadata.position=[2,1];
pard.loadmetadata.TooltipString=sprintf('Load micromanager Metadata.txt file.');
pard.loadmetadata.Optional=true;
pard.calibrate.object=struct('Style','pushbutton','String','auto calibration');
pard.calibrate.position=[2,2];
pard.calibrate.TooltipString=sprintf('calibrate gain and offset from images (from simpleSTORM)');
pard.calibrate.Optional=true;

pard.camparbutton.object=struct('Style','pushbutton','String','set Cam Parameters');
pard.camparbutton.position=[2,3.5];
pard.camparbutton.Width=1.5;
pard.camparbutton.TooltipString=sprintf('Edit camera acquisition parameters.');
pard.plugininfo.type='WorkflowModule'; 
pard.plugininfo.description='Allows editing metadata, or loading from a file.';
end