classdef calibrater3D_so<interfaces.DialogProcessor
    properties
        SXY
    end
    methods
        function obj=calibrater3D_so(varargin)   
            obj@interfaces.DialogProcessor(varargin{:}) ;
%              obj.showresults=true;
%              obj.guiPar.FieldHeight=obj.guiPar.FieldHeight-1;obj.guiPar.Vrim=obj.guiPar.Vrim-20;
        end
        function out=run(obj,p)
            
            fit3ddir=strrep(pwd,'SMAP','fit3D');
            addpath(fit3ddir);
            fit3ddir=strrep(pwd,'SMAP','fit3Dcspline');
            addpath(fit3ddir);
                       
             %get bead positions
             minframes=800/p.dz;
            locDatacopy=obj.locData.copy;
            locDatacopy.regroup(150,minframes/2);
            locg=locDatacopy.getloc({'xnm','ynm','numberInGroup','filenumber'},'layer',1,'Position','roi','removeFilter','filenumber','grouping','grouped');
            beadind=find(locg.numberInGroup>minframes);
            xb=locg.xnm(beadind);
            yb=locg.ynm(beadind);
            filenumber=locg.filenumber(beadind);

            
             locfiles=obj.locData.files.file;

             for k=1:length(locfiles) 
                if isfield(locfiles(k).info,'imagefile')&&~isempty(locfiles(k).info.imagefile)
                    filename=locfiles(k).info.imagefile;
                else
                    filename=locfiles(k).info.basefile;
                end
                il=getimageloader(obj,filename);
                file{k}.name=il.file;
                file{k}.metadata=il.metadata;
                infile=filenumber==k;
                roi=il.metadata.roi;
                pixelsize=il.metadata.cam_pixelsize_um;
                file{k}.posx=xb(infile)/pixelsize(1)/1000-roi(1);
                file{k}.posy=yb(infile)/pixelsize(end)/1000-roi(2);
                il.close;
                posall{k}=horzcat(file{k}.posx,file{k}.posy);
             end
             
             smappos.positions=posall;
             smappos.imageROI=roi;
             cg=calibrate3D_GUI(smappos);
             cg.guihandles.dz.String=num2str(p.dz);
             cg.guihandles.filelist.String=getFieldAsVector(file,'name');
             [p, f]=fileparts(file{1}.name);
             cg.guihandles.outputfile.String=[p fileparts(f) '_3dcal.mat'];
             if length(xb)>0
                 cg.guihandles.posfromsmap.Value=true;
             end
             cg.guihandles.emgain.Value=file{1}.metadata.EMon;
             out=[];
        end
        
%         function initGui(obj)
%             setvisible(0,0,obj)
% %             beaddistribution_callback(0,0,obj)           
%         end
 
        function pard=guidef(obj)
            pard=guidef(obj);
        end
    end
end


function pard=guidef(obj)

pard.dzt.object=struct('String','dz (nm)','Style','text');
pard.dzt.position=[1,1];
pard.dz.object=struct('String','10','Style','edit');
pard.dz.position=[1,2];


pard.inputParameters={'cam_pixelsize_um'};
pard.plugininfo.type='ProcessorPlugin';
pard.plugininfo.description=sprintf('Plugin to calibrate 3D PSF. \n According to Li et al, 2017');
end
