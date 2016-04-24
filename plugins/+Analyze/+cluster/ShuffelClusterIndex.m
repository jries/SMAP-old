classdef ShuffelClusterIndex<interfaces.DialogProcessor
    % ShuffelClusterIndex randomizes the index of clusters. This is useful
    % for color-coded plotting of clusters.
    methods
        function obj=ShuffelClusterIndex(varargin)        
            obj@interfaces.DialogProcessor(varargin{:}) ;
            obj.inputParameters={'sr_pixrec','numberOfLayers','sr_pos','sr_size','layers','sr_layerson'};
            obj.history=true;
            obj.showresults=false;
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
        
        function pardef=guipar(obj)
            pardef.plugininfo.name='Shuffle cluster indices';
            pardef.plugininfo.description='ShuffelClusterIndex randomizes the index of clusters. This is useful for color-coded plotting of clusters.';
            pardef.plugininfo.type='ProcessorPlugin';
        end

    end
end

