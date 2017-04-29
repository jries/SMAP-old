classdef density_calculator<interfaces.DialogProcessor
    % density_calculator looks at the neighborhood and counts number of
    % neighbours. locData.clusterdensity=neighbours
    methods
        function obj=density_calculator(varargin)    
            obj@interfaces.DialogProcessor(varargin{:}) ;
            obj.inputParameters={'sr_pixrec','numberOfLayers','sr_pos','sr_size','layers','sr_layerson'};
            obj.history=true;    
            obj.showresults=false;
        end
        
        function out=run(obj,p)
            out=[];
            obj.locData.sort('filenumber','channel','xnm');
            %problem on filtered data, somehow it gets confused.
            [locs,indin]=obj.locData.getloc({'xnm','ynm','channel','frame','znm','ingrouped','inungrouped'},'position','all');
           %    [locs,indin]=obj.locData.getloc({'xnm','ynm','channel','frame','znm','ingrouped','inungrouped'},'layer',find(p.sr_layerson),'position','all');
           
            x=locs.xnm;
            y=locs.ynm;
            
            dx=p.countingsize_xy;
            dz=p.countingsize_z;
            if ~isempty(locs.znm)
                z=locs.znm;
                if p.countingregion.Value==1 %Gauss
                    neighbours=countneighbours3DGauss(double(x),double(y),double(z),double(dx),double(dz));
                else
                    neighbours=countneighbours3Dcirc(double(x),double(y),double(z),double(dx),double(dz));
                end
            else
                if p.countingregion.Value==1 %Gauss
                    neighbours=countneighbours2DGauss(double(x),double(y),double(dx));
                else
                    neighbours=countneighbours2Dcirc(double(x),double(y),double(dx));
                end
            end
            
            if sum(locs.inungrouped)==length(locs.xnm) %ungrouped data

            else
                neighbours=obj.locData.grouped2ungrouped(locs.ingrouped,neighbours);

            end
            neighbourstot=zeros(length(obj.locData.loc.xnm),1);
            neighbourstot(locs.inungrouped)=neighbours;
            obj.locData.setloc('clusterdensity',single(neighbourstot));
            obj.locData.sort('filenumber','channel','frame');
            obj.locData.regroup;
            obj.setPar('locFields',fieldnames(obj.locData.loc));  
        end
        function pard=guidef(obj)
            pard=guidef;
        end
        

    end
end




function pard=guidef
pard.countingregion.object=struct('String','Gauss wighted counting |circle counting','Style','popupmenu','Value',2);
pard.countingregion.object.TooltipString=sprintf('count Gauss-weighted locs (more accurate) or locs in circle/cylinder (faster)');
pard.countingregion.position=[1,1];
pard.countingregion.Width=2;

pard.texta.object=struct('String','size in x,y (nm)','Style','text');
pard.texta.position=[2,1];
pard.countingsize_xy.object=struct('String','12','Style','edit');
pard.countingsize_xy.position=[2,2];
pard.countingsize_xy.isnumeric=1;
pard.countingsize_xy.object.TooltipString=sprintf('radius of circle or sigma of gauss in lateral direction');


pard.text1.object=struct('String','size in z (nm)','Style','text');
pard.text1.position=[3,1];
pard.countingsize_z.object=struct('String','24','Style','edit');
pard.countingsize_z.position=[3,2];
pard.countingsize_z.isnumeric=1;
pard.countingsize_z.object.TooltipString=sprintf('size of cylinder or sigma of gauss in z direction');
pard.plugininfo.name='density calculator (number of neighbours)';
pard.plugininfo.description= 'density_calculator looks at the neighborhood and counts number of neighbours. If grouped or ungrouped data is used depends on setting in layers.';
pard.plugininfo.type='ProcessorPlugin';
end