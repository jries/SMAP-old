classdef plugin4workflow<interfaces.WorkflowModule
    properties
        subpluginpath
        subplugin
    end
    methods
        function obj=plugin4workflow(varargin)
            obj@interfaces.WorkflowModule(varargin{:})
            obj.inputChannels=1;
        end

        function initGui(obj)
            initGui@interfaces.WorkflowModule(obj);
            thisplugin=plugin(obj.subpluginpath{:});
            thisplugin.handle=obj.handle;
            thisplugin.attachLocData(obj.locData);
            thisplugin.attachPar(obj.P);
            if isa(thisplugin, 'interfaces.DialogProcessor')
                obj.guiPar.Vrim=0;
                thisplugin.processorgui=false;
            end
            thisplugin.setGuiAppearence(obj.guiPar);
            thisplugin.makeGui;
            
            obj.subplugin=thisplugin;

        end
        function prerun(obj,p)
            p=obj.getGuiParameters;
           
        end
        function out=run(obj,data,p)
            out=[];
            if ~data.eof
                return
            end
            if isempty(obj.locData.loc)
                disp('no localization data present')
            end
            module=obj.subplugin;
            if isa(module, 'interfaces.DialogProcessor')
                module.processgo;
            elseif isa(module, 'interfaces.GuiModuleInterface')
                disp('no dialog processor. implement with run')
            end
            obj.output(data)
        end

        

    end
end