classdef correct3Daberrations<interfaces.DialogProcessor
    properties
        SXY
    end
    methods
        function obj=correct3Daberrations(varargin)   
            obj@interfaces.DialogProcessor(varargin{:}) ;
             obj.showresults=true;
             obj.history=true;
%              obj.guiPar.FieldHeight=obj.guiPar.FieldHeight-1;obj.guiPar.Vrim=obj.guiPar.Vrim-20;
        end
        function out=run(obj,p)
             obj.setPar('undoModule','correct3Daberrations');
            notify(obj.P,'backup4undo');
%             p.EMon=obj.locData.files.file(1).info.EMon;
            p.RIM=obj.locData.history{1}.children.fitparamters.fitterGUI.children.MLE_GPU_Yiming.refractive_index_mismatch;
%             p.dz=p.dz*p.RIM;
            Zcorr=load(p.calfile);
            Zcorr=Zcorr.ZcorrInterp;
            z=obj.locData.loc.znm;
            zg=obj.locData.grouploc.znm;
            
            dz=Zcorr.interp(ones(size(z))*p.objectivepos,z/p.RIM)*p.RIM;
            dzg=Zcorr.interp(ones(size(zg))*p.objectivepos,zg/p.RIM)*p.RIM;
            obj.locData.loc.znm=z+dz;
            obj.locData.grouploc.znm=zg+dzg;
           
           
        end
        
        function initGui(obj)
%             setvisible(0,0,obj)
%             beaddistribution_callback(0,0,obj)           
        end
 
        function pard=guidef(obj)
            pard=guidef(obj);
        end
    end
end

function loadcall(a,b,obj)
file=obj.getSingleGuiParameter('calfile');
if isempty(file)
    file=obj.getPar('lastSMLFile');
end
[f,p]=uigetfile(file);
obj.setGuiParameters(struct('calfile',[p f]));
end

function pard=guidef(obj)
pard.objectivepost.object=struct('String','distance objective above glass (nm)','Style','text');
pard.objectivepost.position=[1,1];
pard.objectivepost.Width=2;
pard.objectivepos.object=struct('String','0','Style','edit');
pard.objectivepos.position=[1,3];
pard.objectivepos.Width=1;


pard.calfilet.object=struct('String','calibration file','Style','text');
pard.calfilet.position=[2,1];
pard.calfilet.Width=3;
pard.calfile.object=struct('String','','Style','edit');
pard.calfile.position=[3,1];
pard.calfile.Width=3;

pard.loadbutton.object=struct('Style','pushbutton','String','Load','Callback',{{@loadcall,obj}}); 
pard.loadbutton.position=[3,4];
pard.loadbutton.Width=1;

pard.inputParameters={'cam_pixelsize_nm'};
pard.plugininfo.type='ProcessorPlugin';


end
