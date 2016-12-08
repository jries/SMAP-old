classdef imageloaderMM<interfaces.imageloaderSMAP
    %imageloaderMM image loader for micromanager  tiff stack files
    %   Detailed explanation goes here
    
    properties
        reader
    end
    
    methods
        function obj=imageloaderMM(varargin)
            obj@interfaces.imageloaderSMAP(varargin{:});
        end
        function open(obj,file)
            initMM(obj);
            obj.reader = javaObjectEDT('org.micromanager.acquisition.TaggedImageStorageMultipageTiff',fileparts(file), false, [], false, false, true);
            obj.file=file;
%             obj.reader=bfGetReader(file);
            md=obj.getmetadata;
            [p,f]=fileparts(file);
            obj.metadata.basefile=[p filesep f];
            
        end
        function image=getimage(obj,frame)
            image=readstack(obj,frame);
        end
        
        function allmd=getmetadatatags(obj)
            img=obj.reader;
            imgmetadata=img.getImageTags(0,0,0,0);
            summarymetadata=img.getSummaryMetadata;
            
            allmd=gethashtable(imgmetadata);
            alls=gethashtable(summarymetadata);
            try
            comments=char(img.getDisplayAndComments.get('Comments'));
            allmd(end+1,:)={'Comments direct',comments};
            catch
            end
            %direct
            try
            troi=textscan(imgmetadata.get('ROI'),'%d','delimiter','-');
            allmd(end+1,:)={'ROI direct',num2str(troi{:}')};
            catch err
            end
            framesd=img.lastAcquiredFrame;
            allmd(end+1,:)={'frames direct',num2str(framesd)};
            
            allmd=vertcat(allmd,alls);
            obj.allmetadatatags=allmd;
                
        
        end
        
    end
    
end

function initMM(obj)
% dirs={'ij.jar'
% 'plugins/Micro-Manager/MMAcqEngine.jar'
% 'plugins/Micro-Manager/MMCoreJ.jar'
% 'plugins/Micro-Manager/MMJ_.jar'
% 'plugins/Micro-Manager/clojure.jar'
% 'plugins/Micro-Manager/bsh-2.0b4.jar'
% 'plugins/Micro-Manager/swingx-0.9.5.jar'
% 'plugins/Micro-Manager/swing-layout-1.0.4.jar'
% 'plugins/Micro-Manager/commons-math-2.0.jar'
%  'plugins/Micro-Manager/ome-xml.jar'
%  'plugins/Micro-Manager/scifio.jar'
%  'plugins/Micro-Manager/guava-17.0.jar'
%  'plugins/Micro-Manager/loci-common.jar'
%  'plugins/Micro-Manager/slf4j-api-1.7.1.jar'};

%     if ispc
%         MMpath='C:/Program Files/Fiji/scripts';
%     else
%         MMpath='/Applications/Fiji.app/scripts';
%     end
try
obj.createGlobalSetting('MMpath','Directories2','The directory of Micro-Manager/ in w hich ij.jar is found:',struct('Style','dir','String','MMpath'))
MMpath=obj.getGlobalSetting('MMpath');  
if ~exist(MMpath,'dir')       
    errordlg('cannot find Micro-Manager, please select Micro-Manager directory in menu SMAP/Preferences/Directotries2...')
    return
end


% for k=1:length(dirs)
%     dirs{k}=[MMpath filesep strrep(dirs{k},'/',filesep)];
% end
plugindir=[MMpath filesep 'plugins' filesep 'Micro-Manager' filesep];
allf=dir([plugindir '*.jar']);
dirs={allf(:).name};

for k=1:length(dirs)
    dirs{k}=[MMpath filesep 'plugins' filesep 'Micro-Manager' filesep strrep(dirs{k},'/',filesep)];
end

dirs{end+1}=  [MMpath filesep 'ij.jar']; 
javaaddpath(dirs);
catch err
    disp('java path not added, as imageloader was not called from SMAP. add manually to javaclasspath');
end
end


function image=readstack(obj,imagenumber)
img=obj.reader.getImage(0,0,imagenumber,0);
if isempty(img)
    image=[];
    return
end
image=img.pix;
image=reshape(image,obj.metadata.Width,obj.metadata.Height)';

%    if imagenumber<=obj.reader.getImageCount()
%        image=bfGetPlane(obj.reader,imagenumber);
%    else
%        image=[];
%    end
end

