classdef generalStatistics<interfaces.SEEvaluationProcessor
    methods
        function obj=generalStatistics(varargin)        
                obj@interfaces.SEEvaluationProcessor(varargin{:});
        end
        function out=run(obj,p)
            layerson=p.sr_layerson;%obj.locData.parameters.sr_layerson;
            out.PSFlayers=[];
            out.locplayers=[];
            out.Nlayers=[];
%             roisize=obj.site.sePar.Settings.siteroi/2;
            roisize=p.se_siteroi;
%             if obj.display
%              h=obj.setoutput('p1');
%              h.NextPlot='replace';
%             end
            for k=1:p.numberOfLayers
                if layerson(k)
                    locs=obj.getLocs({'locprecnm','PSFxnm','xnm','ynm'},'layer',k,'size',[roisize,roisize]/2);  
                    
                    if isfield(locs,'PSFxnm')
                    psf=mean(locs.PSFxnm);
                    out.PSFlayers(end+1)=psf;
                    out.(['layers' num2str(k)]).PSF=psf;
                    end
                    
                    locp=mean(locs.locprecnm);
                    out.locplayers(end+1)=locp;
                    out.(['layers' num2str(k)]).locp=locp;
                    
                    N=length(locs.xnm);
                    out.Nlayers(end+1)=N;
                    out.(['layers' num2str(k)]).N=N;
                     
                end
            end
            
            out.Nch=[];
            out.Nchg=[];
            
            maxc=4;
            locs=obj.getLocs({'channel'},'size',[roisize,roisize]/2,'grouping','ungrouped'); 
            glocs=obj.getLocs({'channel'},'size',[roisize,roisize]/2,'grouping','grouped');           
            for c=0:maxc
                ind=locs.channel==c;
                Nch=sum(ind);
                if Nch>0
                    out.Nch(end+1)=Nch;
                    out.(['channels' num2str(c)]).N=Nch;
                    
                    Nchg=sum(glocs.channel==c);
                    out.Nchg(end+1)=Nchg;
                    out.(['channels' num2str(c)]).Ng=Nchg;
                end          
            end   
        end
        function pard=guidef(obj)
            pard=guidef;
        end
    end
end

function pard=guidef
pard.inputParameters={'numberOfLayers','sr_layerson','se_cellfov','se_sitefov','se_siteroi'};
pard.plugininfo.type='ROI_Evaluate';
end