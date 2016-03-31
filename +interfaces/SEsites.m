classdef SEsites<matlab.mixin.Copyable
    properties
        pos
        ID
        info
        annotation
        evaluation
        name
        image
        handles 
        sePar
        indList
    end
    methods
        function obj=SEsites(varargin)
            obj@matlab.mixin.Copyable;
            if nargin>0
                sein=varargin{1};
                fn=properties(sein);
                for k=1:length(fn)
                    obj.(fn{k})=sein.(fn{k});
                end
            end
        
        list.value=1;
        list.string={'empty'};
        obj.annotation.list1=list;
        obj.annotation.list2=list;
        obj.annotation.list3=list;
        obj.annotation.list4=list;
        obj.annotation.comments='';
        
        line.pos=zeros(2);line.value=0;line.angle=0;line.length=0;
        obj.annotation.line1=line;
        obj.annotation.line2=line;
        obj.annotation.rotationpos=line;


        end
    end
end