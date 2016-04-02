classdef Loader_csv<interfaces.DialogProcessor
    methods
        function obj=Loader_csv(varargin)        
                obj@interfaces.DialogProcessor(varargin{:}) ;
                obj.inputParameters={'mainGui'};
        end
        
        function out=load(obj,p,file,mode)
            if nargin<4
                mode=getfilemode(file);
            end
            loadfile(obj,p,file,mode);
        end
        function pard=pardef(obj)
            pard=pardef;
        end
        function run(obj,p)
            [f,p]=uigetfile(obj.info.extensions);
            obj.load(p,[p f]);
            initGuiAfterLoad(obj);
        end
        function clear(obj,file,isadd)
            if isadd 
                obj.locData.clear('filter');
            else
                obj.locData.clear;
            end
        end
    end
end




function pard=pardef
info.name='CSV loader';
info.extensions={'*.csv';'*.*'};
info.dialogtitle='select any CSV file';
pard.plugininfo=info;     
end

function loadfile(obj,p,file,mode)
%look for file description
path=fileparts(file);
descfile=dir([path filesep 'file-description.xml']);
if ~isempty(descfile)
    pfile=pxml2p(xml2struct([path filesep descfile(1).name]));
else
    pfile.frame=1;pfile.xnano=2;pfile.ynano=3;pfile.znano=4;pfile.intensity=5;
end
dat=csvread(file,1,0);
sdat=size(dat);
% filedat=load(file);
% filedat.filename=file;
filenumber=obj.locData.files.filenumberEnd+1;

locData=interfaces.LocalizationData;
locData.addloc('frame',dat(:,pfile.frame));
locData.addloc('xnm',single(dat(:,pfile.xnano)));
locData.addloc('ynm',single(dat(:,pfile.ynano)));
if isfield(pfile,'znano')
locData.addloc('znm',single(dat(:,pfile.znano)));
end
locData.addloc('phot',single(dat(:,pfile.intensity)));

phot=single(dat(:,pfile.intensity));
zd=zeros(sdat(1),1,'single');
locData.addloc('bg',zd);

psfnm=150;psfznm=500;
locData.addloc('locprecznm',psfznm./sqrt(phot));
locData.addloc('locprecnm',psfnm./sqrt(phot));
locData.addloc('PSFxnm',zd+psfnm);
locData.addloc('channel',zd);
locData.addloc('filenumber',zeros(sdat(1),1,'uint8')+filenumber);


obj.locData.addLocData(locData);

filestruc=locData.files.file;
filestruc.name=file;
filestruc.info=struct('Width',256,'Height',256,'roi',[0 0 256 256],'pixsize',psfnm/1000);
if obj.locData.files.filenumberEnd==0
    obj.locData.files.file=filestruc;
else
    obj.locData.files.file(filenumber)=copyfields(obj.locData.files.file(1),filestruc,fieldnames(obj.locData.files.file(1)));
end
obj.locData.files.filenumberEnd=length(obj.locData.files.file);

end

function p=pxml2p(pxml)
fna=fieldnames(pxml.description);
fn=setdiff(fna,{'separator'});
pn=copyfields([],pxml.description,fn);
for k=1:length(fn)
    p.(fn{k})=str2double(pn.(fn{k}).Text)+1;
end
end
