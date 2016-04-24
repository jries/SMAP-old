classdef BlinkingMoviePresentation<interfaces.DialogProcessor
    methods
        function obj=BlinkingMoviePresentation(varargin)        
                obj@interfaces.DialogProcessor(varargin{:});
            obj.inputParameters={'sr_size','sr_pos','sr_pixrec'};
        end
        
        function out=run(obj,p)           
            file=obj.locData.files.file(p.dataselect.Value);
            locs=obj.locData.getloc({'xnm','ynm','frame','locprecnm'},'layer',1,'position','roi');
            makeBlinkMovie(locs,file,p);

        end
        function pard=guidef(obj)
            pard=guidef;
        end
    end
end




function pard=guidef

pard.texta.object=struct('Style','text','String','filter from layer 1');
pard.texta.position=[1,3];

pard.dataselect.object=struct('Style','popupmenu','String','File');
pard.dataselect.position=[1,1];

pard.texta.object=struct('Style','text','String','number of frames');
pard.texta.position=[2,1];


pard.numberOfFrames.object=struct('Style','edit','String','1000');
pard.numberOfFrames.position=[2,2];

pard.textb.object=struct('Style','text','String','first frame');
pard.textb.position=[3,1];

pard.frame_min.object=struct('Style','edit','String','1');
pard.frame_min.position=[3,2];

pard.removedark.object=struct('Style','checkbox','String','remove dark frames','Value',1);
pard.removedark.position=[4,1];
%sum

pard.outputFormat.object=struct('Style','popupmenu','String',{{'MPEG-4','Uncompressed AVI','Motion JPEG 2000'}});
pard.outputFormat.position=[5,1];
pard.outputFormat.Width=2;
pard.syncParameters={{'filelist_short','dataselect',{'String'}}};
pard.plugininfo.name='BlinkingMoviePresentation';

end