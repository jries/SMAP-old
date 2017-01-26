classdef GainOffsetFFT<interfaces.DialogProcessor
    methods
        function obj=GainOffsetFFT(varargin)        
                obj@interfaces.DialogProcessor(varargin{:});
                obj.showresults=true;
        end
        function initGui(obj)
                fn=obj.getPar('loc_filename');
                if ~isempty(fn)
                    obj.guihandles.imagefile.String=fn;
                end
        end
        
        function out=run(obj,p)
            fn=p.imagefile;
            ax=obj.initaxis('image');
            
            il=imageloaderAll(fn);
            if numel(p.numberimages)==1
                p.numberimages=1:p.numberimages;
            end
            img=il.getmanyimages(p.numberimages,'mat');
            imgh=double(img);
            imagesc(imgh(:,:,1));
            RNStd=0;
            ax=obj.initaxis('result');
            [gain,off]=pcfo_j((imgh), p.kthresh, RNStd, 1, 1, p.tiles);
%             fi=obj.getPar('loc_fileinfo');
            fi=il.getmetadata;
            if p.EMon
                fak=p.emgain*2;
            else
                fak=1;
            end
            gainb=fak/gain;
            obj.setGuiParameters(struct('gain',gainb,'offset',off));
            
            disp(['conversion: ', num2str(gainb), ', set: ', num2str(fi.conversion), ', offset: ', num2str(off), ', set: ', num2str(fi.offset)])
            
            out=[];
        end
        function pard=guidef(obj)
            pard=guidef(obj);
        end
        function loadfile(obj,file)
            ht='imagefile';
            il=imageloaderAll(file);
            fi=il.getmetadata;
            obj.setGuiParameters(struct(ht,file,'emgain',fi.emgain,'EMon',double(fi.EMon)))
    
        end
    end
end


function load_callback(a,b,obj)
fn=obj.getSingleGuiParameter('imagefile');
if ~exist(fn,'file')
fn=obj.getPar('loc_filename');
end
if isempty(fn)
    fn='*.tif';
end
[f,p]=uigetfile(fn);
if f
    obj.loadfile(file);
end


end

function setpar_callback(a,b,obj)
p=obj.getGuiParameters;
settings=obj.getPar('loc_fileinfo');
% if p.fixoffsetc
%     p.offset=p.fixoffset;
% end
settings=copyfields(settings,struct('conversion',p.gain,'offset',p.offset,'emgain',p.emgain,'EMon',p.EMon));
obj.setPar('loc_fileinfo',settings);
end

function pard=guidef(obj)

%EMon , emgain: set here. load-button: overwrites with own metadata.

l1=2;
pard.t1.object=struct('Style','text','String','Parameters PCFO:');
pard.t1.position=[l1,1];
pard.t1.Width=1;

pard.kthresht.object=struct('Style','text','String','kthresh');
pard.kthresht.position=[l1,2];
pard.kthresht.Width=.5;
pard.kthresh.object=struct('Style','edit','String','0.9');
pard.kthresh.position=[l1,2.5];
pard.kthresh.Width=.3;

pard.tilest.object=struct('Style','text','String','Tiles');
pard.tilest.position=[l1,3];
pard.tilest.Width=.5;
pard.tiles.object=struct('Style','edit','String','3');
pard.tiles.position=[l1,3.5];
pard.tiles.Width=.3;


pard.numberimagest.object=struct('Style','text','String','# images');
pard.numberimagest.position=[l1,4];
pard.numberimagest.Width=.5;
pard.numberimages.object=struct('Style','edit','String','1');
pard.numberimages.position=[l1,4.5];
pard.numberimages.Width=.5;

pard.imagefilet.object=struct('Style','text','String','Image:');
pard.imagefilet.position=[1,1];
pard.imagefilet.Width=.5;
pard.imagefile.object=struct('Style','edit','String','');
pard.imagefile.position=[1,1.5];
pard.imagefile.Width=3;
pard.imagefileb.object=struct('Style','pushbutton','String','load','Callback',{{@load_callback,obj,1}});
pard.imagefileb.position=[1,4.5];
pard.imagefileb.Width=.5;

% pard.darkfilet.object=struct('Style','checkbox','String','Dark image');
% pard.darkfilet.position=[2,1];
% pard.darkfilet.Width=1;
% pard.darkfile.object=struct('Style','edit','String','');
% pard.darkfile.position=[2,2];
% pard.darkfile.Width=2.5;
% pard.darkfileb.object=struct('Style','pushbutton','String','load','Callback',{{@load_callback,obj,2}});
% pard.darkfileb.position=[2,4.5];
% pard.darkfileb.Width=.5;

% pard.fixoffsetc.object=struct('Style','checkbox','String','Fix offset:');
% pard.fixoffsetc.position=[3,1];
% pard.fixoffsetc.Width=1;
% pard.fixoffset.object=struct('Style','edit','String','0');
% pard.fixoffset.position=[3,2];
% pard.fixoffset.Width=1;

pard.EMon.object=struct('Style','checkbox','String','EM gain:');
pard.EMon.position=[3,3];
pard.EMon.Width=1;
pard.emgain.object=struct('Style','edit','String','100');
pard.emgain.position=[3,4];
pard.emgain.Width=1;

l2=5;

pard.offsett.object=struct('Style','text','String','Offset');
pard.offsett.position=[l2,1];
pard.offsett.Width=.6;
pard.offset.object=struct('Style','edit','String','0');
pard.offset.position=[l2,1.6];
pard.offset.Width=0.4;

pard.gaint.object=struct('Style','text','String','conversion');
pard.gaint.position=[l2,2];
pard.gaint.Width=.6;
pard.gain.object=struct('Style','edit','String','1');
pard.gain.position=[l2,2.6];
pard.gain.Width=0.4;

pard.setcampar.object=struct('Style','pushbutton','String','set Cam Parameter','Callback',{{@setpar_callback,obj}});
pard.setcampar.position=[l2,3];
pard.setcampar.Width=2;

% pard.syncParameters={{'loc_filename','imagefile',{'String'}}};
pard.plugininfo.type='ProcessorPlugin';
end