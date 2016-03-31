classdef TifSaver<interfaces.DialogProcessor
    methods
        function obj=TifSaver(varargin)  
                obj@interfaces.DialogProcessor(varargin{:}) ;
                obj.inputParameters={'filelist_long','mainfile','mainGui','numberOfLayers','sr_layerson'};
        end
        
        function out=save(obj,p)
            obj.status('save tif file')
            fn=p.filelist_long.selection;
            [path,file]=fileparts(fn);
            of=[path filesep file '.tif'];
            
            [f,path]=uiputfile(of);
            if f
                img=obj.getPar('sr_image');
                imwrite(img.image,[path f]);
            end
            obj.status('save done')
          
        end
%         function pard=pardef(obj)
%             pard=pardef;
%         end
        function run(obj,p)
            obj.save(p)
        end        

    end
end