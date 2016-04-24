classdef Locstatistics<interfaces.DialogProcessor
    methods
        function obj=Locstatistics(varargin)           
            obj@interfaces.DialogProcessor(varargin{:}) 
            obj.inputParameters={'sr_layerson'};
            obj.showresults=true;
        end
        
        function out=run(obj,p)
            

        fields={'filenumber','frame','phot','locprecnm','znm','PSFxnm','locprecznm','numberInGroup','bg'};
        if p.useroi
            position='roi';
        else
            position='all';
        end



%         usefields={{'photons','Nloc'},{'photons','mu'},{'lifetime','mu'},{'background','mean'},{'PSFxnm','max'},{'frames','falloff'}};
        layers=find(p.sr_layerson);
            if p.filter
                for m=length(layers):-1:1
                    locs{m}=obj.locData.getloc(fields,'layer',layers(m),'position',position);
                end
            else
                locs{2}=obj.locData.getloc(fields,'position',position,'grouping','grouped');
                locs{1}=obj.locData.getloc(fields,'position',position,'grouping','ungrouped');
            end
            out=make_statistics2(locs,p,true);

%             make_statistics2(obj.locData,p)
%             out=0;
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
pard.checkphot.position=[3,1];
pard.checkphot.Width=2;

pard.photrange.object=struct('String','0','Style','edit');
pard.photrange.position=[3,3];

pard.plugininfo.name='statistics';
end