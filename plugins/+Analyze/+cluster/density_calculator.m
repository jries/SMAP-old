classdef density_calculator<interfaces.DialogProcessor
    methods
        function obj=density_calculator(varargin)    
            obj@interfaces.DialogProcessor(varargin{:}) ;
        end
        
        function out=run(obj,p)
            out=[];
            obj.locData.sort('filenumber','channel','xnm');
            x=obj.locData.loc.xnm;
            y=obj.locData.loc.ynm;
            
            dx=p.countingsize_xy;
            dz=p.countingsize_z;
            if isfield(obj.locData.loc,'znm')
                z=obj.locData.loc.znm;
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
            obj.locData.setloc('neighbours',single(neighbours));
            obj.locData.sort('filenumber','channel','frame');
            obj.locData.regroup;
            obj.setPar('locFields',fieldnames(obj.locData.loc));  
        end
        function pard=pardef(obj)
            pard=pardef;
        end
        

    end
end




function pard=pardef


pard.countingregion.object=struct('String','Gauss wighted counting |circle counting','Style','popupmenu','Value',2);
pard.countingregion.object.TooltipString=sprintf('count Gauss-weighted locs (more accurate) or locs in circle/cylinder (faster)');
pard.countingregion.position=[1,1];



pard.texta.object=struct('String','size in x,y (nm)','Style','text');
pard.texta.position=[2,1];
pard.countingsize_xy.object=struct('String','4','Style','edit');
pard.countingsize_xy.position=[2,2];
pard.countingsize_xy.isnumeric=1;
pard.countingsize_xy.object.TooltipString=sprintf('radius of circle or sigma of gauss in lateral direction');




pard.text1.object=struct('String','size in z (nm)','Style','text');
pard.text1.position=[3,1];
pard.countingsize_z.object=struct('String','8','Style','edit');
pard.countingsize_z.position=[3,2];
pard.countingsize_z.isnumeric=1;
pard.countingsize_z.object.TooltipString=sprintf('size of cylinder or sigma of gauss in z direction');
pard.plugininfo.name='density_calculator';

end