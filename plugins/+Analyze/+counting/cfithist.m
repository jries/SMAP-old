classdef cfithist<interfaces.DialogProcessor
    methods
        function obj=cfithist(varargin)        
                obj@interfaces.DialogProcessor(varargin{:});
        end
        
        function out=run(obj,p)
            histogram=obj.getResults('counting_histogram');
%             histogram=obj.locData.guiData.counting.histogram;
            pout=cluster_mmaple_fithist(p,histogram);
           
%             locs=obj.locData.getloc({'frame','xnm','ynm','phot','bg','PSFxnm','locprecnm'},'layer',1,'position','roi');

%             pout.par=cluster_counting(locs,p);
            obj.setGuiParameters(pout);
            out.histogram=histogram;
            out.fit=pout;
        end
        
        function refit_callback(obj)
        end
        function pard=pardef(obj)
            pard=pardef(obj);
        end
        

    end
end




function pard=pardef(obj)
pard.text1.object=struct('String','parameters','Style','text');
pard.text1.position=[3,1];

pard.N0_fit.object=struct('String','N0','Style','radiobutton');
pard.N0_fit.position=[2,1];

pard.N0_v.object=struct('String','10','Style','edit');
pard.N0_v.position=[2,2];
pard.N0_v.isnumeric=1;


pard.pmature_fit.object=struct('String','p mature','Style','radiobutton');
pard.pmature_fit.position=[3,1];

pard.pmature_v.object=struct('String','.5','Style','edit');
pard.pmature_v.position=[3,2];
pard.pmature_v.isnumeric=1;


pard.pblink_fit.object=struct('String','p blink','Style','radiobutton');
pard.pblink_fit.position=[4,1];

pard.pblink_v.object=struct('String','.2','Style','edit');
pard.pblink_v.position=[4,2];
pard.pblink_v.isnumeric=1;


pard.monomer_fit.object=struct('String','monomer fraction','Style','radiobutton');
pard.monomer_fit.position=[5,1];

pard.monomer_v.object=struct('String','.2','Style','edit');
pard.monomer_v.position=[5,2];
pard.monomer_v.isnumeric=1;




pard.ncluster_fit.object=struct('String','n in cluster','Style','radiobutton');
pard.ncluster_fit.position=[6,1];

pard.ncluster_v.object=struct('String','0','Style','edit');
pard.ncluster_v.position=[6,2];
pard.ncluster_v.isnumeric=1;

pard.blinkmode.object=struct('Style','popupmenu','String','Poisson|Exponential|Noblink');
pard.blinkmode.position=[7,1];
pard.blinkmode.Width=2;


pard.t1.object=struct('String','fitrange','Style','text');
pard.t1.position=[2,3];

pard.fitrange_min.object=struct('String','4','Style','edit');
pard.fitrange_min.position=[2,4];
pard.fitrange_min.Width=0.5;
pard.fitrange_max.object=struct('String','60','Style','edit');
pard.fitrange_max.position=[2,4.5];
pard.fitrange_max.Width=0.5;

pard.t2.object=struct('String','fit:','Style','text');
pard.t2.position=[4,3];
pard.t2.Width=0.25;
pard.fitselection.object=struct('String',{{'cumulative distribution','histogram','weighted histogram'}},'Style','popupmenu');
pard.fitselection.position=[4,3.25];
pard.fitselection.Width=1.75;

pard.plugininfo.name='fit brightness histogram';

% pard.text2.object=struct('String','length scale nm','Style','text');
% pard.text2.position=[8,1];
% 
% pard.lengthscale.object=struct('String','15','Style','edit');
% pard.lengthscale.position=[8,2];
% pard.lengthscale.isnumeric=1;

% pard.results1='results 1';
% pard.results3='other results';


end