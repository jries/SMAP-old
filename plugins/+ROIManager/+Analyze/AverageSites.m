classdef AverageSites<interfaces.DialogProcessor&interfaces.SEProcessor
    methods
        function obj=AverageSites(varargin)        
                obj@interfaces.DialogProcessor(varargin{:});
            obj.inputParameters={'se_viewer'};
        end
        
        function out=run(obj,p)  
            out=[];
            
            obj.locData.addfile(p.name);
            initGuiAfterLoad(obj);
            obj.SE.processors.preview.updateFilelist;
            newfile=obj.locData.files.filenumberEnd;
            locnew=obj.locData.loc;
            sites=obj.locData.SE.sites;
            used=false(size(locnew.xnm));
            for k=1:length(sites)
                [locs,indloc]=obj.locData.getloc('xnm','position',sites(k));
                locnew.xnm(indloc)=locnew.xnm(indloc)-sites(k).pos(1);
                locnew.ynm(indloc)=locnew.ynm(indloc)-sites(k).pos(2);
                locnew.filenumber(indloc)=newfile;
                used=used|indloc;
            end
            fn=fieldnames(locnew);
            for k=1:length(fn)
                obj.locData.addloc(fn{k},locnew.(fn{k})(used))
            end
            obj.locData.regroup;
            obj.locData.filter;
            
            
            %try: add empty file, there put averaged sites
            %for every site: loc.xnm-site.pos(1)+xpossite
           
        end
        function pard=guidef(obj)
            pard=guidef;
        end
    end
end




function pard=guidef

pard.t1.object=struct('String','average sites','Style','text');
pard.t1.position=[1,1];
pard.t1.Width=1;

pard.post.object=struct('String','position (x,y)','Style','text');
pard.post.position=[2,1];
pard.post.Width=1;

pard.pos.object=struct('String','0,0','Style','edit');
pard.pos.position=[2,2];
pard.pos.Width=1;


pard.sortselection.object=struct('String',{{'all','use'}},'Style','popupmenu');
pard.sortselection.position=[2,3];
pard.sortselection.Width=1;

pard.namet.object=struct('String','name','Style','text');
pard.namet.position=[3,1];
pard.namet.Width=1;

pard.name.object=struct('String','average','Style','edit');
pard.name.position=[3,2];
pard.name.Width=1;


pard.plugininfo.type='ROI_Analyze';


end