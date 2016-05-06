classdef AnalyzeRingCME<interfaces.DialogProcessor&interfaces.SEProcessor
    properties
        sites
        results
        resultsfigure
    end
    methods
        function obj=AnalyzeRingCME(varargin)        
            obj@interfaces.DialogProcessor(varargin{:});
            obj.inputParameters={};
        end
        
        function out=run(obj,p)  
            if p.restrictcheck
                if length(p.usesites)==1
                    range=1:min(p.usesites,length(obj.SE.sites));
                else
                    range=p.usesites;
                end
            else 
                range=1:length(obj.SE.sites);
            end
            obj.sites=obj.SE.sites(range);
                
            if isempty(obj.sites)
                disp('no sites loaded')
                return
            end
            obj.results=cmeresults(obj.sites);   
            obj.resultsfigure=cmeRingplotresults(p,obj.results);
            out=obj.results;

        end
        function pard=guidef(obj)
            pard=guidef(obj);
        end

        function save_callback(obj,a,b)
            p=obj.getAllParameters;
            selection=p.savemenu.selection;
%             [~,~,ext]=fileparts(selection);
            [file,path]=uiputfile(selection);
            if path
                switch selection
                    case 'results.mat'
                        results=obj.results;
                        save([path file],'results');
                    case 'sites.mat'
                        sites=obj.sites;
                        save([path file],'sites');
                    case 'results.pdf'
                        export_fig([path file],'-pdf','-nocrop',obj.resultsfigure)
%                         saveas(obj.resultsfigure,[path file],'epsc')

                    case 'average.tif'
                        image=obj.results.sumimage/obj.results.numsites;
                        saveastiff(uint16(image/max(image(:))*2^16),[path file])
                    case 'rad_av.mat'
                        results.sumimage=obj.results.sumimage;
                        results.sumrdensity=obj.results.sumrdensity;
                        results.numberofsites=obj.results.numsites;
                        save([path file],'results');
                    case 'all (not sites)'
                        results=obj.results;
                        save([path 'results.mat'],'results');
                        
                        results.sumimage=obj.results.sumimage;
                        results.sumrdensity=obj.results.sumrdensity;
                        results.numberofsites=obj.results.numsites;
                        save([path 'rad_av.mat'],'results');
                        
                        image=obj.results.sumimage/obj.results.numsites;
                        saveastiff(uint16(image/max(image(:))*2^16),[path file])
                         export_fig([path file],'-pdf','-nocrop',obj.resultsfigure)

                end
            end
        end
    end
end

function results=cmeresults(sites)
circfitfields={'evaluation','CME2DRing','circfit'};
imfitfields={'evaluation','CME2DRing','imfit'};
results.N=getFieldAsVector(sites,circfitfields{:},'Ncirc1');
results.rc=getFieldAsVector(sites,circfitfields{:},'r1');

results.images=getFieldAsVector(sites,imfitfields{:},'image');
results.dr=getFieldAsVector(sites,imfitfields{:},'dr1');
results.ro=getFieldAsVector(sites,imfitfields{:},'r1');

results.sigma=getFieldAsVector(sites,imfitfields{:},'sigma1');

results.rdensity=getFieldAsVector(sites,imfitfields{:},'profiles1','rdensity');
results.ac=getFieldAsVector(sites,imfitfields{:},'profiles1','thetaAC');
results.rdensityn=getFieldAsVector(sites,imfitfields{:},'profiles1','rn');
results.acthetan=getFieldAsVector(sites,imfitfields{:},'profiles1','thetan');

results.filenumber=getFieldAsVector(sites,'info','filenumber');
results.sitenames=getFieldAsVector(sites,'name');
    
    
filenumberrange=1:max(results.filenumber);
results.filenumberrange=filenumberrange;
Nfilemean=0;
Nfilemedian=0;
% Nnorm=N;
N=results.N;
Nmed=median(N);
Nmean=mean(N);
for k=filenumberrange
    ink=results.filenumber==k;
    Nfilemean(k)=mean(N(ink));
    Nfilemedian(k)=median(N(ink));
    Nnormmean(ink)=N(ink)/Nfilemean(k)*Nmean;
    Nnormmedian(ink)=N(ink)/Nfilemedian(k)*Nmed;  
end

results.Nfilemean=Nfilemean;
results.Nfilemedian=Nfilemedian;
results.Nnormmean=Nnormmean;
results.Nnormmedian=Nnormmedian;

results.numsites=length(sites);
sumim=zeros(size(results.images{1}));
sumrdensity=zeros(size(results.rdensity{1}));
for k=1:length(sites)
%     size(results.images{k})
    sumim=results.images{k}+sumim;
    sumrdensity=results.rdensity{k}+sumrdensity;
end
results.sumimage=sumim;
results.sumrdensity=sumrdensity;

results.sumrdensityn=results.rdensityn{1};

end


function pard=guidef(obj)


pard.savebutton.object=struct('String','Save','Style','pushbutton','Callback',{{@obj.save_callback}});
pard.savebutton.position=[1,1];

pard.savemenu.object=struct('String',{{'all (not sites)','results.pdf','average.tif','results.mat','sites.mat','rad_av.mat'}},'Style','popupmenu');
pard.savemenu.position=[2,1];

pard.t1.object=struct('String','max dr/ro','Style','text');
pard.t1.position=[1,3];
pard.maxdrro.object=struct('String','1.5','Style','edit');
pard.maxdrro.position=[1,4];

pard.restrictcheck.object=struct('String','only sites','Style','checkbox','Value',0);
pard.restrictcheck.position=[2,3];
pard.usesites.object=struct('String','150','Style','edit');
pard.usesites.position=[2,4];
pard.usesites.TooltipString='first N sites, or LIst of sites or s1:s2 notation';
pard.plugininfo.type='ROI_Analyze';
end