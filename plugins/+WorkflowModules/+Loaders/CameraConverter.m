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
        offset
        adu2phot
    end
    methods
        function obj=CameraConverter(varargin)
            obj@interfaces.WorkflowModule(varargin{:})
            obj.loc_cameraSettings=interfaces.metadataSMAP;
            fn=fieldnames(obj.loc_cameraSettingsStructure);
            for k=1:length(fn)
                obj.loc_cameraSettings.assigned.(fn{k})=false;
            end
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
        function setmetadata(obj,overwrite)
            if nargin<2
                overwrite=false;
            end
            md=obj.loc_cameraSettings;
            obj.loc_cameraSettings=interfaces.metadataSMAP;
            obj.loc_cameraSettings=copyfields(obj.loc_cameraSettings,md);
            if ~overwrite
                settings=obj.getPar('loc_fileinfo');
                fn=fieldnames(settings);
                for k=1:length(fn)
                    if settings.assigned.(fn{k})
                        obj.loc_cameraSettings.(fn{k})=settings.(fn{k});
                    end
                end
            end
            obj.setPar('loc_cameraSettings',obj.loc_cameraSettings);
            obj.setPar('EMon',obj.loc_cameraSettings.EMon);
            obj.updatefileinfo;
%             if obj.loc_cameraSettings.EMon
%             	obj.EMexcessNoise=2;
%             else
%                 obj.EMexcessNoise=1;
%             end
        end
        function setcamerasettings(obj,fi)
            fn=fieldnames(fi);
            fn2=properties(obj.loc_cameraSettings);
            fna=intersect(fn,fn2);
            for k=1:length(fna)
                obj.loc_cameraSettings.(fna{k})=fi.(fna{k});
            end
        end
        function updatefileinfo(obj)
           fi=obj.getPar('loc_fileinfo');
           md2=obj.loc_cameraSettings;
           fn=fieldnames(fi);fn2=properties(md2);
           fna=intersect(fn,fn2);
           fi=copyfields(fi,md2,fna);
           obj.setPar('loc_fileinfo',fi);
        end
        function prerun(obj,p)
            pc=obj.loc_cameraSettings;
            if ~pc.EMon
                pc.emgain=1;
            end
            obj.offset=pc.offset;
            obj.adu2phot=(pc.conversion/pc.emgain);
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
           
            imphot=single((data.data-obj.offset))*obj.adu2phot;
               datao=data;
               datao.data=imphot;  
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
fi=obj.getPar('loc_fileinfo');
for k=length(fn):-1:1
    fields{k}=fn{k};
    if isfield(fi,fn{k})
        defAns{k}=num2str(fi.(fn{k}));
    else
        defAns{k}=num2str(obj.loc_cameraSettings.(fn{k}));
    end
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
obj.setmetadata(true);
%     if obj.loc_cameraSettings.EMon
%         obj.EMexcessNoise=2;
%     else
%         obj.EMexcessNoise=1;
%     end
%     obj.globpar.parameters.loc_cameraSettings=obj.loc_cameraSettings; %doesnt work
end
end

function loadmetadata_callback(a,b,obj)

% return
finf=obj.getPar('loc_fileinfo');
if ~isempty(finf)
    ft=[(finf.basefile) filesep '*.*'];
else
    ft='*.*';
end

[f,p]=uigetfile(ft,'Select metadata.txt, tif or _sml.mat');
[~,~,ext ]=fileparts(f);
switch ext
    case '.mat'
        l=load([p f]);
        metadata=l.saveloc.file(1).info;



    case '.txt'
        par.metadatafile=[p f];
        obj.setGuiParameters(par);
        obj.updateGui;       
        metadata=getmetadataMMtxt([p f]); 
    otherwise
        imloader=imageloaderAll([p f],finf);
        metadata=imloader.metadata;
        metadata.allmetadata=metadata;
end
        obj.setcamerasettings(metadata);
        obj.updatefileinfo;
% obj.readmetadata;
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

pard.text.object=struct('Style','text','String','Acquisition parameters:');
pard.text.position=[1,1];
pard.text.Width=1.5;
pard.text.Optional=true;
% pard.metadatafile.object=struct('Style','edit','String',' ');
% pard.metadatafile.position=[2,1];
% pard.metadatafile.Width=4;
% pard.metadatafile.Optional=true;
pard.loadmetadata.object=struct('Style','pushbutton','String','Load metadata');
pard.loadmetadata.position=[1,2.5];
pard.loadmetadata.TooltipString=sprintf('Load micromanager Metadata.txt file.');
pard.loadmetadata.Optional=true;
% pard.calibrate.object=struct('Style','pushbutton','String','auto calibration');
% pard.calibrate.position=[2,2];
% pard.calibrate.TooltipString=sprintf('calibrate gain and offset from images (from simpleSTORM)');
% pard.calibrate.Optional=true;

pard.camparbutton.object=struct('Style','pushbutton','String','set Cam Parameters');
pard.camparbutton.position=[1,3.5];
pard.camparbutton.Width=1.5;
pard.camparbutton.TooltipString=sprintf('Edit camera acquisition parameters.');
pard.plugininfo.type='WorkflowModule'; 
pard.plugininfo.description='Allows editing metadata, or loading from a file.';
end