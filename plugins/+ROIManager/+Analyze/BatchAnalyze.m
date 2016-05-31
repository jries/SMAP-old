classdef BatchAnalyze<interfaces.DialogProcessor&interfaces.SEProcessor
    properties
    end
    methods
        function obj=BatchAnalyze(varargin)        
            obj@interfaces.DialogProcessor(varargin{:});
            obj.inputParameters={};
        end
        function initGui(obj)
%             g=obj.getPar('mainGui');
            plugins=obj.getPar('menu_plugins');
            fn=fieldnames(plugins.ROIManager.Analyze);
            obj.guihandles.analyzer.String=fn;
        end
        
        function out=run(obj,p)  
            out=[];
            g=obj.getPar('mainGui');
            
            [f,pfad]=uigetfile('*.mat');
            gf=g.children.guiFile;
            gf.loadbutton_callback(0,0,0,pfad,f);
            
            if p.redrawall
                g.locData.SE.processors.preview.redrawall;
            end
            analzyers=fieldnames(g.children.guiSites.children.Analyze.children);
            selectedan=p.analyzer.selection;
            if any(strcmp(analzyers,selectedan))
                analyzering=g.children.guiSites.children.Analyze.children.(selectedan);
                analyzering.processgo;
            else
                disp('selected analyzer should be added to ROImanager/Analyzers')
            end
            
            
            g.locData.savelocs([pfad 'rois_' f]);

            [~,filen]=fileparts([pfad f]);
            analyzering.saveall(pfad,filen);
        end
        function pard=guidef(obj)
            pard=guidef(obj);
        end
    end
end


function pard=guidef(obj)


pard.t1.object=struct('String','select analyzer','Style','text');
pard.t1.position=[1,1];
pard.analyzer.object=struct('String','empty','Style','popupmenu');
pard.analyzer.position=[1,2];
pard.analyzer.Width=3;

pard.redrawall.object=struct('String','redraw all','Style','checkbox','Value',1);
pard.redrawall.position=[2,1];

pard.plugininfo.type='ROI_Analyze';
end