classdef showdriftinfo<interfaces.DialogProcessor
    methods
        function obj=showdriftinfo(varargin)        
                obj@interfaces.DialogProcessor(varargin{:}) ;
               
        end
        
        function out=run(obj,p)
            out=[];
            driftinfo=obj.locData.files.file(p.dataselect.Value).driftinfo;
            dx=driftinfo.dx;
            dy=driftinfo.dy;
            dxplot=driftinfo.dxplot;
            dyplot=driftinfo.dyplot;
            dxt=driftinfo.dxt;
            dyt=driftinfo.dyt;
            binframes=driftinfo.binframes;
            framesall=1:length(dxt);
            results_ax2=initaxis(p.resultstabgroup,'drift vs frame');
            subplot(1,2,1)
            hold off
            plot(dxplot)
            hold on
            plot(dx,'k','LineWidth',1.5);
            sx=(max(dx)-min(dx));
            ylim([min(dx)-sx/2 max(dx)+sx/2])
            axis tight

            subplot(1,2,2)
            hold off
            plot(dyplot)
            hold on
            plot(dy,'k','LineWidth',1.5);

            sy=(max(dx)-min(dx));
            ylim([min(dy)-sy/2 max(dy)+sy/2])
            axis tight
            
            initaxis(p.resultstabgroup,'drift vs frame  final');

            hold off
            plot(binframes,dx,'x',framesall,dxt,'k')
            hold on
            plot(binframes,dy,'o',framesall,dyt,'r')
            xlabel('frame')
            ylabel('dx, dy (nm)')
            drawnow

            initaxis(p.resultstabgroup,'dx vs dy');
            hold off
            plot(dxt,dyt,'k')
            hold on
            plot(dx,dy,'ro')
            plot(dx(1),dy(1),'gx')
            xlabel('dx')
            ylabel('dy')
            drawnow
            axis equal
          
        end
        function pard=pardef(obj)
            pard=pardef;
        end
    end
end




function pard=pardef
pard.dataselect.object=struct('Style','popupmenu','String','File');
pard.dataselect.position=[2,1];
pard.syncParameters={{'filelist_short','dataselect',{'String'}}};
end