classdef simpleParameterInterface<handle
    properties
        inputParameters
        parameters
    end
    methods
        function setParameters(obj,varargin)
            for k=1:length(varargin)
                obj.parameters=copyfields(obj.parameters,varargin{k},obj.inputParameters);
            end
        end
        function p=getAllParameters(obj)
            p=obj.parameters;
        end
    end
end
