classdef Loader_sml<interfaces.DialogProcessor;
    methods
        function obj=Loader_sml(varargin)        
                obj@interfaces.DialogProcessor(varargin{:}) ;
                obj.inputParameters={'mainGui'};
        end
        
        function out=load(obj,p,file,mode)
            if nargin<4
                mode=getfilemode(file);
            end
            if isempty(p)
                p.updateGuiPar=false;
            end
            loadfile(obj,p,file,mode);
        end
        function run(obj,p)
            [f,p]=uigetfile(obj.info.extensions);
            obj.load(p,[p f]);
            initGuiAfterLoad(obj);
        end
        function pard=guidef(obj)
            pard=guidef;
        end
        function clear(obj,file,isadd)
            if nargin<3
                isadd=false;
            end
            if isadd 
                obj.locData.clear('filter');
            else
                obj.locData.clear;
            end
        end        
    end
end




function pard=guidef
info.name='SML loader';
info.extensions={'*.mat';'*.*'};
info.dialogtitle='select any SMLM file (_sml or _fitpos)';
pard.plugininfo=info;
pard.plugininfo.type='LoaderPlugin';
end

function loadfile(obj,p,file,mode)
% fobj=obj.locData.files;
% obj.locData.files.filenumberEnd=obj.locData.files.filenumberEnd+1; %write back to .files
filedat=load(file);
filedat.filename=file;

filenumber=obj.locData.files.filenumberEnd;
switch mode
    case 'sml'
        [templocData,GUIsettings,siteexplorer]=load_smlV3(filedat);
        obj.setPar('lastSMLFile',file);
    case 'fitpos'
        templocData=loadfitposV2(filedat);
        GUIsettings=[];
        siteexplorer=[];
    case 'sites'
        [templocData,siteexplorer]=load_sites(filedat);
        GUIsettings=[];
    otherwise 
        disp('file format not recognized');
        return;
end
indout=templocData.loc.locprecnm>250|imag(templocData.loc.locprecnm)~=0|isnan(templocData.loc.locprecnm);
templocData.removelocs(indout);

%correct filenumber: .loc, files.filenumber, add files.
fn=fieldnames(templocData.loc);
ldat=length(templocData.loc.(fn{1}));
if ~isfield(templocData.loc,'filenumber')
    templocData.loc.filenumber=ones(ldat,1);
end
if ~isfield(templocData.loc,'channel')
    templocData.loc.filenumber=ones(ldat,1,'single');
end
nfiles=length(templocData.files.file);
templocData.loc.filenumber=templocData.loc.filenumber+filenumber;
obj.locData.addLocData(templocData);
his=templocData.history;
obj.locData.history(end+1:end+length(his))=his;
for k=1:nfiles
%     templocData.files.file(k).number=templocData.files.file(k).number+filenumber;
    templocData.files.file(k).number=k+filenumber;
    if ~myisfield(templocData.files.file(k).info,'roi')||isempty(templocData.files.file(k).info.roi)||~isnumeric(templocData.files.file(k).info.roi)
        templocData.files.file(k).info.roi=[0 0 templocData.files.file(k).info.Width templocData.files.file(k).info.Height];
    end
    if ~isfield(templocData.files.file(k).info,'cam_pixelsize_um')
        templocData.files.file(k).info.cam_pixelsize_um=templocData.files.file(k).info.pixsize([1 1]);
    end
end

if isfield(templocData.files.file(k).info,'pixsize')
    templocData.files.file(k).info=rmfield(templocData.files.file(k).info,'pixsize');
end

newfilenumbers=filenumber+1:filenumber+nfiles;
filestruc=templocData.files.file;
if obj.locData.files.filenumberEnd==0
    obj.locData.files.file=filestruc;
else
    for k=1:length(filestruc)    

        obj.locData.files.file(filenumber+k)=copyfields(obj.locData.files.file(1),filestruc(k),fieldnames(obj.locData.files.file(1)));
    end
end
obj.locData.files.filenumberEnd=length(obj.locData.files.file);

if ~isempty(GUIsettings) %write back parameters
%     button=questdlg('Restore saved GUI Parameters?','GUI parameters','Yes','No','No');
    if p.updateGuiPar %strcmpi(button,'Yes')
        if isfield(GUIsettings,'par')
            GUIsettings=convertparameters(GUIsettings);
        end
        p.mainGui.setGuiParameters(GUIsettings,true)
    end
end

if isempty(siteexplorer)||siteexplorer.numberOfFiles==0
    siteexplorer=interfaces.SiteExplorer;
    for k=1:length(templocData.files.file)
        siteexplorer.addFile(templocData.files.file(k).name,templocData.files.file(k).number,templocData.files.file(k).info)
    end
%     newfilenumbers=1:templocData.files.file(end).number;
else
    rgp=obj.getPar('ROI_restorparamters');
    if ~p.updateGuiPar && length(siteexplorer.sites)>0 && (isempty(rgp) ||rgp)%restore parameters only if not already globally restored
        psites=GUIsettings.children.guiSites;
        p.mainGui.children.guiSites.setGuiParameters(psites,true);
    end
end
se=obj.locData.SE;
   se.addSites(siteexplorer,newfilenumbers, templocData.files.file)
   
   
end
