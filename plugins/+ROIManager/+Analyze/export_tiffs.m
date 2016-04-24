classdef export_tiffs<interfaces.DialogProcessor&interfaces.SEProcessor
    methods
        function obj=export_tiffs(varargin)        
                obj@interfaces.DialogProcessor(varargin{:});
            obj.inputParameters={};
        end
        
        function out=run(obj,p)  
            sites=obj.SE.sites;
            if p.export_selected
                sev=obj.getPar('se_viewer');
                pv=sev.getAllParameters;
                selectedsites=pv.sitelist.Value;
            else
              
                selectedsites=1:length(sites);
            end
            mainfile=obj.getPar('mainfile');
            path=fileparts(mainfile);
            prefix='img_.tif';
            [f,path]=uiputfile([path filesep prefix]);
            
            if ~f
                return
            end
            [~,f]=fileparts(f);
     
            for k=selectedsites
                site=sites(k);
                imold=site.image.image;
                site.image=[];
                site.image=obj.SE.plotsite(site,-1);
                site.image.image=imold;
                options.color=true;
                options.comp='lzw';
                filen=[f strrep(site.name,'.','_')];
                fhere= [filen '.tif'];
                imout=uint8(site.image.image*255);
                
                saveastiff(imout,[path fhere],options);
                if length(site.image.layers)>1
                    for ll=1:length(site.image.layers)
                        
                        iml=site.image.layers(ll).images.renderimages.image;
                        filell= [filen '_' int2str(ll) '.tif'];
                        imoutll=uint8(iml*255);
                        saveastiff(imoutll,[path filell],options);
                    end
                end
                site.image.layers=[];site.image.composite=[];
            end
            out=0;
        end
        function pard=guidef(obj)
            pard=guidef;
        end
    end
end




function pard=guidef

pard.export_selected.object=struct('String','exprot only selected sites','Style','checkbox');
pard.export_selected.position=[1,1];
pard.export_selected.Width=2;

end