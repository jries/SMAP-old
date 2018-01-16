classdef export_3Dvolumes<interfaces.DialogProcessor&interfaces.SEProcessor
    methods
        function obj=export_3Dvolumes(varargin)        
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
            
            if p.export_cells
                cells=obj.SE.cells;
                for k=1:length(cells)
                    cell=cells(k);
                    imold=cell.image.image;
                    cell.image=[];
                    cell.image=obj.SE.plotsite(cell,-1);
                    cell.image.image=imold;
                    options.color=true;
                    options.comp='lzw';
                    filen=[f 'cell_C' num2str(cell.ID)  '_F' num2str(cell.info.filenumber)];
                    fhere= [filen '.tif'];
                    imout=uint8(cell.image.image*255);

                    saveastiff(imout,[path fhere],options);
                end
            end
            
            out=0;
        end
        function pard=guidef(obj)
            pard=guidef;
        end
    end
end




function pard=guidef

pard.export_selected.object=struct('String','export only selected sites','Style','checkbox');
pard.export_selected.position=[1,1];
pard.export_selected.Width=2;

pard.pixrect.object=struct('String','pixelsize (nm) [x y z]','Style','text');
pard.pixrect.position=[2,1];
pard.pixrect.Width=1;

pard.pixrec.object=struct('String','7','Style','edit');
pard.pixrec.position=[2,2];
pard.pixrec.Width=1;

pard.sizerect.object=struct('String','size (pixels) [x y z]','Style','text');
pard.sizerect.position=[3,1];
pard.sizerect.Width=1;

pard.sizerec.object=struct('String','50','Style','edit');
pard.sizerec.position=[3,2];
pard.sizerec.Width=1;

pard.plugininfo.type='ROI_Analyze';

end