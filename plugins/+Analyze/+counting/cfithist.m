classdef cfithist<interfaces.DialogProcessor
    %  Copyright (c)2017 Ries Lab, European Molecular Biology Laboratory,
    %  Heidelberg. This file is part of Single Molecule Analysis Platform (SMAP).
    
    % cfithist fits a histogram 
    
    methods
        function obj=cfithist(varargin)        
                obj@interfaces.DialogProcessor(varargin{:});
        end
        
        function out=run(obj,p)
            
            if p.bootstrap
                numberOfSubsets = length(obj.getResults('counting_histogram'));
            else
                numberOfSubsets = 1;
            end
            
            mature = [];
            blink = [];
            for i = 1:numberOfSubsets
                if p.bootstrap
                    allHistogram=obj.getResults('counting_histogram');
                    histogram=allHistogram{i};
                else
                    histogram=obj.getResults('counting_histogram');
                end
%               histogram=obj.locData.guiData.counting.histogram;
                pout_ref=cluster_mmaple_fithist(p,histogram);
                
                while 1
                    pout=cluster_mmaple_fithist(pout_ref,histogram);

                    if round(pout.pmature_v,3) == round(pout_ref.pmature_v,3) & round(pout.pblink_v,3) == round(pout_ref.pblink_v,3)
                        break
                    else
                        pout_ref=pout;
                    end
                end
                blink(i) = pout_ref.pblink_v;
                mature(i) = pout_ref.pmature_v;
                
%               locs=obj.locData.getloc({'frame','xnm','ynm','phot','bg','PSFxnm','locprecnm'},'layer',1,'position','roi');

%               pout.par=cluster_counting(locs,p);
            end
            
            if p.bootstrap
                pout.blinkStd = std(blink); pout.matureStd = std(mature);
                pout.blinkMean = mean(blink); pout.matureMean = mean(mature);
            end
            obj.setGuiParameters(pout);
            out.histogram=histogram;
            out.fit=pout;
        end
        
        function refit_callback(obj)
        end
        function pard=guidef(obj)
            pard=guidef(obj);
        end
        

    end
end




function pard=guidef(obj)
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

pard.bootstrap.object=struct('String','Bootstrap','Style','checkbox', 'Value', 0);
pard.bootstrap.position=[5,3];
pard.bootstrap.Width=1;

pard.t_bsMean.object=struct('String','Mean','Style','text');
pard.t_bsMean.position=[6,4];
pard.t_bsMean.Width=0.5;

pard.t_bsStd.object=struct('String','Std','Style','text');
pard.t_bsStd.position=[6,4.5];
pard.t_bsStd.Width=0.5;

pard.t_bsMature.object=struct('String','Mature','Style','text');
pard.t_bsMature.position=[7,3.5];
pard.t_bsMature.Width=0.5;

pard.t_bsBlink.object=struct('String','Blink','Style','text');
pard.t_bsBlink.position=[8,3.5];
pard.t_bsBlink.Width=0.5;

pard.matureMean.object=struct('String','-','Style','edit');
pard.matureMean.position=[7,4];
pard.matureMean.Width=0.5;

pard.blinkMean.object=struct('String','-','Style','edit');
pard.blinkMean.position=[8,4];
pard.blinkMean.Width=0.5;

pard.matureStd.object=struct('String','-','Style','edit');
pard.matureStd.position=[7,4.5];
pard.matureStd.Width=0.5;

pard.blinkStd.object=struct('String','-','Style','edit');
pard.blinkStd.position=[8,4.5];
pard.blinkStd.Width=0.5;


pard.plugininfo.name='fit brightness histogram';
pard.plugininfo.type='ProcessorPlugin';
% pard.text2.object=struct('String','length scale nm','Style','text');
% pard.text2.position=[8,1];
% 
% pard.lengthscale.object=struct('String','15','Style','edit');
% pard.lengthscale.position=[8,2];
% pard.lengthscale.isnumeric=1;

% pard.results1='results 1';
% pard.results3='other results';


end