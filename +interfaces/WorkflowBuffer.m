classdef WorkflowBuffer<handle
    properties(Access=private)
        dat={};
        currenttag=1;
        bufferincrease=10;
        complete=false;
    end
    methods
        function add(obj,data,tag)
            lod=length(obj.dat);
            if lod<tag
                    obj.dat{tag}=data;
                obj.complete(lod+1:tag+obj.bufferincrease)=false;
            else
            
            ld=length(obj.dat{tag});
 
            obj.dat{tag}(ld+1)=data;
            end
            if data.numberInTag>0 && length(obj.dat{tag})>=data.numberInTag
                obj.complete(tag)=true;
            end
            if tag>obj.currenttag
                if data.numberInTag==0&&length(obj.dat{obj.currenttag})>=1
                    obj.complete(obj.currenttag)=true;
                end
                obj.currenttag=tag;    
            end
            if data.eof
                 obj.complete(tag)=true;
            end
  
        end
        function out=iscomplete(obj,tag)
            if nargin<2
                out=obj.complete;
            else
                out=obj.complete(tag);
            end
        end
        function data=get(obj,tag)
            data=obj.dat{tag}(1:end);
%             obj.dat{tag}.empty;
%              obj.dat{tag}=interfaces.WorkflowData;
            obj.dat{tag}.data=[];
            obj.complete(tag)=false;          
        end
    end
end