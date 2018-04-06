classdef gate
    properties
        name;
        indices;
        channelNames;
        fullFilename; %opt
        community_labels %opt
    end
    
methods
    
    function obj = gate(name, indices, channelNames, ...
            fullFilename, community_labels)
        obj.name = name;
        obj.indices = indices;
        obj.channelNames = channelNames;
        if nargin > 3
            obj.fullFilename = fullFilename;
        end
        
        if nargin > 4
            obj.community_labels = community_labels;
        end

    end
end
end
        