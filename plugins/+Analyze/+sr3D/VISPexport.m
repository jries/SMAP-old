classdef VISPexport<interfaces.DialogProcessor

    methods
        function obj=VISPexport(varargin)        
            obj@interfaces.DialogProcessor(varargin{:}) ;
             obj.inputParameters={'sr_layerson'};
        end
        function pard=guidef(obj)
            pard.plugininfo.type='ProcessorPlugin';
        end
        function out=run(obj,p)
            locs=obj.locData.getloc({'xnm','ynm','znm','locprecznm','locprecnm','phot','frame'},'position','roi');
            [path,f]=fileparts(obj.locData.files.file(1).name);
            [f,path]=uiputfile([path filesep f '.3dlp']);
            [~,fx]=fileparts(f);
            facfwhm=2.3;
            for k=1:5
                if p.sr_layerson(k)
                    fout=[path filesep fx '_l' num2str(k) '.3dlp'];
                    outmatrix=[locs.xnm-mean(locs.xnm) locs.ynm-mean(locs.ynm) locs.znm locs.locprecnm*facfwhm locs.locprecnm*facfwhm locs.locprecznm*facfwhm   locs.phot locs.frame];
                    dlmwrite(fout,outmatrix,'delimiter','\t')
                end
            end          
            out=0;
        end
    end
end