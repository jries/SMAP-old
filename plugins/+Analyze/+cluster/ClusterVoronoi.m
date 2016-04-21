classdef ClusterVoronoi<interfaces.DialogProcessor
    methods
        function obj=ClusterVoronoi(varargin)        
            obj@interfaces.DialogProcessor(varargin{:}) ;
            obj.inputParameters={'sr_pixrec','numberOfLayers','sr_pos','sr_size','layers','sr_layerson'};
            obj.history=true;
        end
        function out=run(obj,p)
            
            out=vrender(obj,p);
        end
        function out=renderfunction(obj,locs,p)
             pos=[p.sr_pos(1)-p.sr_size(1) p.sr_pos(2)-p.sr_size(2) 2*p.sr_size(1) 2*p.sr_size(2)];
                out=vrender(locs,p.sr_pixrec,pos);
        end
        function pard=pardef(obj)
            pard=pardef;
        end
    end
end
function out=vrender(obj,p)
out=[];
                warning('off','MATLAB:scatteredInterpolant:DupPtsAvValuesWarnId')
                locs=obj.locData.getloc({'xnm','ynm','phot','frame','znm','ingrouped','inungrouped'},'layers',find(p.sr_layerson),'position','roi');
                A=zeros(length(locs.xnm),9);

                A(:,2)=locs.frame;
                A(:,4)=locs.xnm-min(locs.xnm);
                A(:,5)=locs.ynm-min(locs.ynm);  
                A(:,7)=locs.phot;
                if ~isempty(locs.znm)
                    A(:,6)=locs.znm;
                end
                Vor = VorArea( A, 1, size(A, 1) );
                Area = Vor{1};
                ic = Vor{5};
                Areacorr = zeros(size(A,1),1);
                Areacorr(:) = Area(ic,:);

%                 
                if sum(locs.inungrouped)==length(locs.xnm) %ungrouped data
%                     
                ico=zeros(length(locs.inungrouped),1,'single');
                ico(locs.inungrouped)=ic;
                obj.locData.setloc('clusterindex',ico);
                
                ceo=zeros(length(locs.inungrouped),1,'single');
                ceo(locs.inungrouped)=Areacorr;
                obj.locData.setloc('clusterdensity',ceo);
%                 end
                obj.setPar('locFields',fieldnames(obj.locData));
                obj.locData.regroup;
                else %if clustering performed on grouped data -> write cluster data to all ungrouped data 
                    %sort according to groupindex
                    
                end

                    
                

end
function pard=pardef
pard.t1.object=struct('String','Voronoi cluster analysis from SharpViSu','Style','text');
pard.t1.position=[1,2];
pard.t1.Width=4;


end

