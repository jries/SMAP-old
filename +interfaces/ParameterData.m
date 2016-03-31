classdef ParameterData<handle
    %global parameter object shared by all modules
    %due to performance: use of global variable. This means: only one copy
    %of SMAP can be open at a single time.
    properties
        par %dummy, taken out for performance and replaced by global variable
        globalSettings
        globalSettingsFile='settings/temp/globalsettings.txt';
        results
        autoresults
    end
    events
        sr_render %tells GuiRender to render again: update the image
        sr_display %tells GuiRender to display again (e.g. layers switched off) without re-rendering
        loc_initialize %initializes workflow files. 
        %TODO take out completely. 
        backup4undo %backs up current obj.locData 
    end
    methods
        function set.par(obj,v) %for performance. take out st.par, get.par in the future when handle classes are faster
            % called when P.par=v is called
            global SMAPparameters
            SMAPparameters=copyfields(SMAPparameters,v);
        end
        function v=get.par(obj)
            % called when v=P.par is called
            global SMAPparameters
            v=SMAPparameters;
        end
        function clear(obj)
            %clears object
            global SMAPparameters
            SMAPparameters=[];
        end
        function loadGlobalSettings(obj)
            if exist(obj.globalSettingsFile,'file')
            obj.globalSettings=readstruct(obj.globalSettingsFile,{},false);
            end
        end
        

    end
    
end