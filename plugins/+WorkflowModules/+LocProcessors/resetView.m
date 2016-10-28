classdef resetView<interfaces.WorkflowModule;
    properties

        
    end
    methods
       function obj=resetView(varargin)
            obj@interfaces.WorkflowModule(varargin{:})
            obj.inputChannels=1; 
        end
        function pard=guidef(obj)
            pard=guidef;
        end
        function initGui(obj)
            initGui@interfaces.WorkflowModule(obj);
        end
        function prerun(obj,p)
           
        end
        function output=run(obj,data,p)
             si=obj.getPar('sr_sizeRecPix');
              if ~isempty(obj.locData.loc)
                mx=myquantilefast(obj.locData.loc.xnm,[0.9995,0.0005],100000);
                maxx=mx(1);minx=mx(2);
                my=myquantilefast(obj.locData.loc.ynm,[0.9995,0.0005],100000);
                maxy=my(1);miny=my(2);
              else
                  disp('cannot find size of image, no reset')
                  return
              end
              obj.setPar('sr_pos',[(maxx+minx)/2 (maxy+miny)/2]);
              pixrec=round(max((maxx-minx)/si(1),(maxy-miny)/si(2)));
              obj.setPar('sr_pixrec',pixrec);
              output=data;
                    
        end
    end
end


function pard=guidef

pard.plugininfo.type='WorkflowModule'; 
pard.plugininfo.description='resets view';
end