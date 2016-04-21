classdef ShuffelClusterIndex<interfaces.DialogProcessor
    methods
        function obj=ShuffelClusterIndex(varargin)        
            obj@interfaces.DialogProcessor(varargin{:}) ;
            obj.inputParameters={'sr_pixrec','numberOfLayers','sr_pos','sr_size','layers','sr_layerson'};
            obj.history=true;
        end
        function out=run(obj,p)
            out=[];
            clusterindex=obj.locData.loc.clusterindex;
            mi=max(clusterindex);
            n=rand(mi,1);
            [~,indsort]=sort(n);
            indsort=vertcat(1,indsort);
            clusterindexnew = indsort(clusterindex+1)-1;
            obj.locData.loc.clusterindex=clusterindexnew;
            obj.locData.regroup;
        end

    end
end

