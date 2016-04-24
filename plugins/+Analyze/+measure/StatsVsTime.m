classdef StatsVsTime<interfaces.DialogProcessor
    methods
        function obj=StatsVsTime(varargin)           
            obj@interfaces.DialogProcessor(varargin{:}) 
            obj.inputParameters={'sr_layerson'};
            obj.showresults=true;
        end
        
        function out=run(obj,p)
            switch p.timefield.selection
                case 'files'
                    fields={'filenumber','frame','phot','locprecnm','znm','PSFxnm','locprecznm','numberInGroup','bg'};
                    if p.useroi
                        position='roi';
                    else
                        position='all';
                    end
                    
                    
                    
                    usefields={{'photons','Nloc'},{'photons','mu'},{'lifetime','mu'},{'background','mean'},{'PSFxnm','max'},{'frames','falloff'},{'locprec','max'},{'locprec','rising'}};
                    layers=find(p.sr_layerson);
                    
                    maxfn=max(obj.locData.getloc('filenumber').filenumber);
                    for filen=1:maxfn
                        if p.filter
                            for m=length(layers):-1:1
                                locs{m}=obj.locData.getloc(fields,'layer',layers(m),'position',position,'filenumber',filen);
                            end
                        else
                            locs{2}=obj.locData.getloc(fields,'position',position,'grouping','grouped','filenumber',filen);
                            locs{1}=obj.locData.getloc(fields,'position',position,'grouping','ungrouped','filenumber',filen);
                        end
                        
                        stats=make_statistics2(locs,p,false);
                        for f=1:length(usefields)
                            fh=usefields{f};
                            stath=stats.(fh{1}).(fh{2});
                            for l=1:length(stath)
                                ps.(fh{1}).(fh{2}){l}(filen)=stath(l);
                            end
                        end
                        
                    end
                    
                    %plot
                    filns=1:maxfn;
                    
                    axall=initaxis(obj.resultstabgroup,'all');
                    delete(axall.Children);
                    hold on
                   ax0=initaxis(obj.resultstabgroup,'falloff');
                    hold off
                    plot(ax0,filns,ps.frames.falloff{1})
                    plot(axall,filns,ps.frames.falloff{1}/max(ps.frames.falloff{1}))
                    
                    for k=2:length(locs)
                    hold on
                    plot(ax0,filns,ps.frames.falloff{k})
                    plot(axall,filns,ps.frames.falloff{k}/max(ps.frames.falloff{k}))
                    end
                    xlabel('filenumber')
                    ylabel('falloff frame')
                    
                    ax1=initaxis(obj.resultstabgroup,'number of localizations');
                    hold off
                    plot(ax1,filns,ps.photons.Nloc{1})
                    plot(axall,filns,ps.photons.Nloc{1}/max(ps.photons.Nloc{1}))
                    for k=2:length(locs)
                    hold on
                    plot(ax1,filns,ps.photons.Nloc{k})
                    plot(axall,filns,ps.photons.Nloc{k}/max(ps.photons.Nloc{k}))
                    end
                    xlabel('filenumber')
                    ylabel('number of localizations')
                    
                    ax2=initaxis(obj.resultstabgroup,'photons (mu)');
                    hold off
                    plot(ax2,filns,ps.photons.mu{1})
                    plot(axall,filns,ps.photons.mu{1}/max(ps.photons.mu{1}))
                    for k=2:length(locs)
                    hold on
                    plot(ax2,filns,ps.photons.mu{k})
                    plot(axall,filns,ps.photons.mu{k}/max(ps.photons.mu{k}))
                    end
                    xlabel('filenumber')
                    ylabel('photons (mu)')
                    
                    ax2b=initaxis(obj.resultstabgroup,'locprec');
                    hold off
                    plot(ax2b,filns,ps.locprec.max{1})
                    hold on
                    plot(ax2b,filns,ps.locprec.rising{1})
                    
                    plot(axall,filns,ps.locprec.rising{1}/max(ps.locprec.rising{1}))
                    plot(axall,filns,ps.locprec.max{1}/max(ps.locprec.max{1}))
                    for k=2:length(locs)
                    hold on
                    plot(ax2b,filns,ps.locprec.max{k})
                    plot(ax2b,filns,ps.locprec.rising{k})
                    plot(axall,filns,ps.locprec.rising{k}/max(ps.locprec.rising{k}))
                    plot(axall,filns,ps.locprec.max{k}/max(ps.locprec.max{k}))                   
                    
                    end
                    xlabel('filenumber')
                    ylabel('photons (mu)')
                    
                    ax3=initaxis(obj.resultstabgroup,'lifetime');
                    hold off
                    plot(ax3,filns,ps.lifetime.mu{1})
                    plot(axall,filns,ps.lifetime.mu{1}/max(ps.lifetime.mu{1}))
                    for k=2:length(locs)
                    hold on
                    plot(ax3,filns,ps.lifetime.mu{k})
                    plot(axall,filns,ps.lifetime.mu{k}/max(ps.lifetime.mu{k}))
                    end
                    xlabel('filenumber')
                    ylabel('lifetime (mu)')
                    
                     ax4=initaxis(obj.resultstabgroup,'BG');
                    hold off
                    plot(ax4,filns,ps.background.mean{1})
                    plot(axall,filns,ps.background.mean{1}/max(ps.background.mean{1}))
                    for k=2:length(locs)
                    hold on
                    plot(ax4,filns,ps.background.mean{k})
                    plot(axall,filns,ps.background.mean{k}/max(ps.background.mean{k}))
                    end
                    xlabel('filenumber')
                    ylabel('mean background')
                    
                    if ~isempty(locs{1}.PSFxnm)
                        ax4=initaxis(obj.resultstabgroup,'PSF');
                        hold off
                        plot(ax4,filns,ps.PSFxnm.max{1})
                        for k=2:length(locs)
                        hold on
                        plot(ax4,filns,ps.PSFxnm.max{k})
                        end
                        xlabel('filenumber')
                        ylabel('max PSFx (nm)')
                    end
                    
                    
                case 'frames'
            end
            out=0;
        end
        function pard=guidef(obj)
            pard=guidef;
        end
    end
end




function pard=guidef
pard.useroi.object=struct('String','use Roi','Style','checkbox','Value',1);
pard.useroi.position=[1,1];

pard.filter.object=struct('String','use layers/filters','Style','checkbox','Value',1);
pard.filter.position=[1,2];

pard.overview.object=struct('String','plot overview','Style','checkbox','Value',0);
pard.overview.position=[1,3];

pard.checkphot.object=struct('String','use manual photon range','Style','checkbox','Value',0);
pard.checkphot.position=[2,1];
pard.checkphot.Width=2;

pard.photrange.object=struct('String','0','Style','edit');
pard.photrange.position=[2,3];

pard.timefield.object=struct('String',{{'files','frames'}},'Style','popupmenu');
pard.timefield.position=[4,1];

pard.t1.object=struct('String','Frame time windows','Style','text');
pard.t1.position=[5,1];
pard.framewindows.object=struct('String','10','Style','edit');
pard.framewindows.position=[5,2];
pard.plugininfo.name='statistics vs frame/file';
end