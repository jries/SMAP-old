classdef PSFmodel<interfaces.GuiModuleInterface
    %UNTITLED2 Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        modelpar
        X
        Y
    end
    
    methods (Abstract=true)
       imout=PSF(obj,x,y,z)

    end
    
end

