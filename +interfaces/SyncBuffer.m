classdef SyncBuffer<handle
    properties(Access=private)
        inputbuffers
        inputchannels
    end
    methods
        function obj=SyncBuffer(inputchannels)
            obj.inputchannels=inputchannels;
            for k=1:inputchannels
                obj.inputbuffers{k}=interfaces.WorkflowBuffer;
            end
        end
        function add(obj,data,tag,inputchannel)
            
            obj.inputbuffers{inputchannel}.add(data,tag);
        end
        function out=iscomplete(obj,tag)
            if nargin<2 %get index for all complete
                
                comp=obj.inputbuffers{1}.iscomplete;
                for k=2:obj.inputchannels
                    comp2=obj.inputbuffers{k}.iscomplete;
                    comp=intersect(comp,comp2);
%                     lm=min(length(comp),length(comp2));
%                     c1=comp(1:lm);c2=comp2(1:lm);
%                     comp=c1&c2;
                end
%                 out=find(comp);
                  out=comp;
            else          
                out=true;
                for k=1:obj.inputchannels
                    out=out&obj.inputbuffers{k}.iscomplete(tag);
                end  
            end
        end
        function out=get(obj,tag)
            for k=1:obj.inputchannels
                out{k}=obj.inputbuffers{k}.get(tag);         
            end
        end
    end
end