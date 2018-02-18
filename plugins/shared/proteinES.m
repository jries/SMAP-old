classdef proteinES

    properties
        %define class properties if needed
        proteinID
        proteinName
        par
    end
    methods
        function obj=proteinES(proteinID, proteinName)   %replace by filename        
            obj.proteinID = proteinID;
            obj.proteinName = proteinName;
            par = table([],[],[],[],[],[],[],[],'VariableNames',{'time', 'dz', 'th', 'de', 'ym', 'shi', 'nm', 'sha'});
            obj.par = par;
        end
        function getTimeImg(t)
            
        end
    end
end


