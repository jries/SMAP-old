classdef AnalyzeNPCquantification<interfaces.DialogProcessor&interfaces.SEProcessor
    properties
    end
    methods
        function obj=AnalyzeNPCquantification(varargin)        
            obj@interfaces.DialogProcessor(varargin{:});
            obj.inputParameters={};
            obj.showresults=true;
        end
        
        function out=run(obj,p)  
            out=runintern(obj,p);
            
        end
        function pard=guidef(obj)
            pard=guidef(obj);
        end

    end
end

function pf=fithist(hi,n,p)

shi=sum(hi);
% hi=hi/shi;

corners=p.corners;rings=p.rings;
n1=find(n>=p.fitrange(1),1,'first');
n2=find(n<=p.fitrange(2),1,'last');
range=n1:n2;
% x=(0:corners)';
x=n(range)';
% clusterfromlabeling(x,corners,rings,.5)
ft=fittype('a*clusterfromlabeling(x,corners,rings,p)','problem',{'corners','rings'});
f=fit(x,hi(range)',ft,'problem',{corners, rings},'Lower',[0 0.01],'Upper',[inf .99],'Start',[shi .1]);
plot(n,f(n),'-g')
plot(x,f(x),'-*r')
pf=f.p;
end

function out=runintern(obj,p)
se=obj.SE;
fields={'evaluation','NPCLabelingQuantify'};
fields2={'evaluation','generalStatistics'};
numfoundint=getFieldAsVector(se.sites,fields{:},'numfoundint');
numfoundrat=getFieldAsVector(se.sites,fields{:},'numfoundrat');
numbercornerassined=getFieldAsVector(se.sites,fields{:},'numbercornerassined');

psf=getFieldAsVector(se.sites,fields2{:},'PSFlayers');
ax0=obj.initaxis('Summary');
axpsf=obj.initaxis('PSF');
histogram(axpsf,psf);
title(axpsf,['PSF range: ' num2str(p.PSFrange)])
xlabel('average PSF (nm)')
indgood=psf>=p.PSFrange(1)&psf<=p.PSFrange(2);
nb=0:p.corners;

numfoundint=numfoundint(indgood);
numfoundrat=numfoundrat(indgood);
numbercornerassined=numbercornerassined(indgood);

ax1=obj.initaxis('from gap');
ax2=axes(ax1.Parent);
subplot(1,2,1,ax1);
hi=hist(numfoundint,nb);
    bar(nb,hi)
    hold on
    pf=fithist(hi,nb,p);
    title(['gap integer: ' num2str(pf,2)])
    axis tight
    results.gapinteger=pf;   
subplot(1,2,2,ax2);
    hold off
    hr=hist(numfoundrat,nb);
    bar(nb,hr)
    hold on
    pf=fithist(hr,nb,p);
    title(['gap fractional: ' num2str(pf,2)])
    axis tight
    results.gapfractional=pf; 

    
 ax3=obj.initaxis('assigned + all');
ax4=axes(ax3.Parent);
subplot(1,2,1,ax3);
    ha=hist(numbercornerassined,nb);
    bar(nb,ha)
    hold on
    pf=fithist(ha,nb,p);
    title(['assigned: ' num2str(pf,2)])
   axis tight
   results.assigned=pf; 
    
subplot(1,2,2,ax4);
    hall=hi+hr+ha;
    bar(nb,hall)
    hold on
    pf=fithist(hall,nb,p);
    title(['all: ' num2str(pf,2)])
    axis tight 
    results.all=pf; 
    

    
    
if isfield(se.sites(1).evaluation.NPCLabelingQuantify,'numcornersfiltered') %not from simulation
    numcorners=getFieldAsVector(se.sites,fields{:},'numcornersfiltered');
    numcornersunf=getFieldAsVector(se.sites,fields{:},'numcorners');
    
    
ax6=obj.initaxis('ground truth');
    hold off
    hnc=hist(numcornersunf,nb);
    bar(nb,hnc)
    hold on
    pf=fithist(hnc,nb);
    results.groundtruth=pf; 
%     
%     hh=clusterfromlabeling(nb,8,2,.4)*sum(hnc);
%     plot(nb,hh)
    title(['true: ' num2str(pf,2)])
    axis tight
ax7=obj.initaxis('comparison');
    plot(numfoundint-numcorners)
    hold on
    plot(numfoundrat-numcorners)
    plot(numbercornerassined-numcorners)
    plot(numcorners/8);
    
    hold off

    legend('integer','prob','assigned','number of corners')

    numfoundintm=mean(numfoundint-numcorners);
    numfoundratm=mean(numfoundrat-numcorners);
    numbercornerassinedm=mean(numbercornerassined-numcorners);

    title(['int: ' num2str(numfoundintm),', rat: ' num2str(numfoundratm), ' ,asssigned: ' num2str(numbercornerassinedm)])
end


dat=struct2table(results);
axp=ax0.Parent;
delete(ax0);
ht=uitable('Parent',axp);
struct2uitable(ht, results,'flip')

out=[];

if p.copy2page
    f=figure;
    ht2=ht.copy;
    ht2.Parent=f;
    axtemp= subplot(3,3,[1 2]);
    ht2.Units='normalized';
    ht2.Position=axtemp.Position;
    delete(axtemp)
%     subplot(2,3,1,axpsf2)
    axpsf2=axpsf.copy;
    axpsf2.Parent=f;
    subplot(3,3,3,axpsf2)
    axis tight
    axt=ax1.copy;
    axt.Parent=f;
    subplot(3,3,4,axt)
    axt=ax2.copy;
    axt.Parent=f;
    subplot(3,3,5,axt)
        axt=ax3.copy;
    axt.Parent=f;
    subplot(3,3,6,axt)
        axt=ax4.copy;
    axt.Parent=f;
    subplot(3,3,7,axt)
    
    if exist(ax6,'var')
        axt=ax6.copy;
    axt.Parent=f;
    subplot(3,3,8,axt) 
    end
end

end

function pard=guidef(obj)
pard.t1.object=struct('String','Corners:','Style','text');
pard.t1.position=[1,1];

pard.corners.object=struct('String','8','Style','edit');
pard.corners.position=[1,2];
pard.corners.Width=0.5;

pard.t2.object=struct('String','Proteins/Corner','Style','text');
pard.t2.position=[2,1];

pard.rings.object=struct('String','4','Style','edit');
pard.rings.position=[2,2];
pard.rings.Width=0.5;


pard.t3.object=struct('String','PSF range','Style','text');
pard.t3.position=[3,1];

pard.PSFrange.object=struct('String','80 150','Style','edit');
pard.PSFrange.position=[3,2];
pard.PSFrange.Width=1;

pard.t4.object=struct('String','fit range histogram','Style','text');
pard.t4.position=[4,1];

pard.fitrange.object=struct('String','3 8','Style','edit');
pard.fitrange.position=[4,2];
pard.fitrange.Width=1;
pard.plugininfo.type='ROI_Analyze';
  

pard.copy2page.object=struct('String','Copy to own page','Style','checkbox');
pard.copy2page.position=[5,1];
pard.copy2page.Width=2;

pard.plugininfo.type='ROI_Analyze';

end