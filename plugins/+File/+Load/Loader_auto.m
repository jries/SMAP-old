classdef Loader_auto<interfaces.DialogProcessor
    methods
        function obj=Loader_auto(varargin)        
                obj@interfaces.DialogProcessor(varargin{:}) ;
                obj.inputParameters={'mainGui'};
        end
        
        function out=load(obj,p,file,mode)
            if nargin<4
                mode=getfilemode(file);
            end
            loadfile(obj,p,file,mode);
        end
        function pard=guidef(obj)
            pard=guidef;
        end
        function run(obj,p)
            [f,p]=uigetfile(obj.info.extensions);
            obj.load(p,[p f]);
            initGuiAfterLoad(obj);
        end
        function clear(obj,file,isadd)
            mode=getfilemode(file);
            loader=getloader(obj,mode);
            loader.clear(file,isadd);
        end
    end  
end




function pard=guidef
info.name='auto loader';
info.extensions={'*.mat;*.tif;*.csv';'*.mat';'*.tif';'*.csv';'*.*'};
info.dialogtitle='select any SMLM position, Tif, csv, settings or workflow file';
pard.plugininfo=info;
end

function loadfile(obj,p,file,mode)            
        loader=getloader(obj,mode);
        loader.load(p,file);
end

        
function loader=getloader(obj,mode)   
        switch mode
            case 'tif'
                loadername='Loader_tif';

            case {'sml','fitpos','sites'}
                loadername='Loader_sml';
                
            case 'guiparameters'
                loadername='Loader_settings';
       
            case 'workflow'
                loadername='Loader_workflow';
            case 'csv'
                loadername='Loader_csv';
              
            otherwise
                disp('file type not recognized')
        end
        loader=plugin('File','Load',loadername,[],obj.P);
        loader.attachLocData(obj.locData);
end
