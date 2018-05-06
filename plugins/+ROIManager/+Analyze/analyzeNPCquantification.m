classdef analyzeNPCquantification<interfaces.DialogProcessor&interfaces.SEProcessor
    properties
    end
    methods
        function obj=analyzeNPCquantification(varargin)        
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

% function pf=fithist(hi,n,p)
% 
% shi=sum(hi);
% % hi=hi/shi;
% 
% corners=p.corners;rings=p.rings;
% n1=find(n>=p.fitrange(1),1,'first');
% n2=find(n<=p.fitrange(2),1,'last');
% range=n1:n2;
% % x=(0:corners)';
% x=n(range)';
% % clusterfromlabeling(x,corners,rings,.5)
% ft=fittype('a*clusterfromlabeling(x,corners,rings,p)','problem',{'corners','rings'});
% f=fit(x,hi(range)',ft,'problem',{corners, rings},'Lower',[0 0.01],'Upper',[inf .99],'Start',[shi .4]);
% plot(n,f(n),'-g')
% plot(x,f(x),'-*r')
% pf=f.p;
% end

function out=runintern(obj,p)
se=obj.SE;
fields={'evaluation','NPCLabelingQuantify'};
fields2={'evaluation','generalStatistics'};

numfoundint=getFieldAsVector(se.sites,fields{:},'numfoundint');
numfoundrat=getFieldAsVector(se.sites,fields{:},'numfoundrat');
numbercornerassinged=getFieldAsVector(se.sites,fields{:},'numbercornerassigned');
numbercornerassigneddirect=getFieldAsVector(se.sites,fields{:},'numbercornerassigneddirect');
filenumber=getFieldAsVector(se.sites,'info','filenumber');
psf=getFieldAsVector(se.sites,fields2{:},'PSFlayers');

use=getFieldAsVector(se.sites,'annotation','use');
numfoundint=numfoundint(use);
numfoundrat=numfoundrat(use);
numbercornerassinged=numbercornerassinged(use);
psf=psf(use);


ax0=obj.initaxis('Summary');
axpsf=obj.initaxis('PSF');
histogram(axpsf,psf);
title(axpsf,['PSF range: ' num2str(p.PSFrange)])
xlabel('average PSF (nm)')
if p.psfcheck
indgood=psf>=p.PSFrange(1)&psf<=p.PSFrange(2);
else
    indgood=true(size(psf));
end

if p.filecheck
indf=false(size(filenumber));
for k=1:length(p.filenumbers)
    indf=indf | filenumber==p.filenumbers(k);
end
indgood=indgood&indf;
end

nb=0:p.corners;

numfoundint=numfoundint(indgood);
numfoundrat=numfoundrat(indgood);
numbercornerassinged=numbercornerassinged(indgood);
numbercornerassigneddirect=numbercornerassigneddirect(indgood);

ax1=obj.initaxis('from gap');
ax2=axes(ax1.Parent);
subplot(1,2,1,ax1);
    p.ploton=false;
    bs_numfoundint=bootstrp(20,@fitNPClabeling,numfoundint,p);
    berr_numfoundint=std(bs_numfoundint);
    p.ploton=true;    
hi=hist(numfoundint,nb);
    bar(nb,hi)
    hold on
    pf=fitNPClabeling(hi,p);
    title(['gap integer: ' num2str(100*pf,'%2.1f') '\pm' num2str(100*berr_numfoundint,'%2.1f')])
    axis tight
    results.gapinteger=pf;   
    results.gapinteger_se=berr_numfoundint;
subplot(1,2,2,ax2);
    p.ploton=false;
    bs_numfoundrat=bootstrp(20,@fitNPClabeling,numfoundrat,p);
    berr_numfoundrat=std(bs_numfoundrat);
    p.ploton=true;    
    
    hold off
    hr=hist(numfoundrat,nb);
    bar(nb,hr)
    hold on
    pf=fitNPClabeling(hr,p);
    title(['gap fractional: ' num2str(100*pf,'%2.1f')  '\pm' num2str(100*berr_numfoundrat,'%2.1f')])
    axis tight
    results.gapfractional=[pf ]; 
    results.gapfractional_se=berr_numfoundrat;
    
 ax3=obj.initaxis('assigned + all');
%  ax3b=axes(ax3.Parent);
ax4=axes(ax3.Parent);
% subplot(1,3,1,ax3);
%     ha=hist(numbercornerassinged,nb);
%     bar(nb,ha)
%     hold on
%     pf=fitNPClabeling(ha,p);
%     title(['assigned: ' num2str(pf,2)])
%    axis tight
%    results.assigned=pf; 

   subplot(1,2,1,ax3);
    ha=hist(numbercornerassigneddirect,nb);
    p.ploton=false;
    bs_assigned=bootstrp(20,@fitNPClabeling,numbercornerassigneddirect,p);
    berr_assigned=std(bs_assigned);
%     [ncbc,bstat]=bootci(30,@fitNPClabeling,numbercornerassigneddirect,p);
%     ncbc=bootci(250,@fitNPClabeling,numbercornerassigneddirect,p);
    
    bar(nb,ha)
    hold on
    p.ploton=true;
    pf=fitNPClabeling(ha,p);
    title(['assigned direct: ' num2str(100*pf,'%2.1f') '\pm' num2str(100*berr_assigned,'%2.1f')])
   axis tight
   results.assigned=pf; 
   results.assigned_se=berr_assigned;
   
subplot(1,2,2,ax4);
    hall=hi+hr+ha;
    bar(nb,hall)
    hold on
    pf=fitNPClabeling(hall,p);
    title(['all: ' num2str(100*pf,'%2.1f')])
    axis tight 
    results.all=pf; 
    

    
    
if isfield(se.sites(1).evaluation.NPCLabelingQuantify,'numcornersfiltered_gt') %not from simulation
    numcorners=getFieldAsVector(se.sites,fields{:},'numcornersfiltered_gt');
    numcornersa=getFieldAsVector(se.sites,fields{:},'numcornersall_gt');
    
    
ax6=obj.initaxis('ground truth');
    p.ploton=false;
    bs_gtf=bootstrp(20,@fitNPClabeling,numcorners(indgood),p);
    berr_gtf=std(bs_gtf);
    p.ploton=true;    
    
    hold off
    hnc=hist(numcorners(indgood),nb);
    bar(nb,hnc)
    hold on
    pf=fitNPClabeling(hnc,p);
    results.groundtruth=pf;
    results.groundtruth_se=berr_gtf;
    title(['true: ' num2str(100*pf,'%2.1f') '\pm' num2str(100*berr_gtf,'%2.1f')])
    axis tight

    ax6b=obj.initaxis('ground truth all');
        p.ploton=false;
    bs_gtfa=bootstrp(20,@fitNPClabeling,numcornersa(indgood),p);
    berr_gtfa=std(bs_gtfa);
    p.ploton=true;  
    
    hold off
    hnc=hist(numcornersa(indgood),nb);
    bar(nb,hnc)
    hold on
    pf=fitNPClabeling(hnc,p);
    results.groundtruthall=pf; 
    results.groundtruthall_se=berr_gtfa;
    title(['true: ' num2str(100*pf,'%2.1f') '\pm' num2str(100*berr_gtfa,'%2.1f')])
    axis tight
    
    rd=(0:length(numcorners)-1)/length(numcorners)/2;
ax7=obj.initaxis('comparison');
    plot(numcorners+rd,numfoundint-numcorners,'+')
    hold on
    plot(numcorners+rd,numfoundrat-numcorners,'x')
%     plot(numcorners+rd,numbercornerassinged-numcorners,'o')
    plot(numcorners+rd,numbercornerassigneddirect -numcorners,'d')
%     plot(numcorners+rd,numbercornerassigneddirect -numcorners,'d')
%     plot(numcorners-p.corners);
    
    hold off

    legend('integer','prob','assigned d','all')

    numfoundintm=mean(numfoundint-numcorners);
    numfoundratm=mean(numfoundrat-numcorners);
    numbercornerassinedm=mean(numbercornerassinged-numcorners);
numbercornerassineddm=mean(numbercornerassigneddirect-numcorners);
    title(['int: ' num2str(numfoundintm,2),', rat: ' num2str(numfoundratm,2), ' ,asssigned d: ' num2str(numbercornerassineddm,2)])
    
    
% ax8=obj.initaxis('comparison h');
% 
% range=-8:8;
% rangep=range(1:end-1)-0.5;
% h1=histcounts(numfoundint-numcorners,range);
% h2=histcounts(numfoundrat-numcorners,range);
% h3= histcounts(numbercornerassinged-numcorners,range);
% stairs(ax8,rangep,h1)
% hold on;
% stairs(rangep,h2);
% stairs(rangep,h3);
% hold off
%   
% 
% axis tight
% xlim(ax8,[-8 8])
% t8={['I:' num2str(mean(numfoundint-numcorners),2) '\pm' num2str(std(numfoundint-numcorners),2)],[ 'p:' num2str(mean(numfoundrat-numcorners),2) '\pm' num2str(std(numfoundrat-numcorners),2)],[ 'a:' num2str(mean(numbercornerassinged-numcorners),2) '\pm' num2str(std(numbercornerassinged-numcorners),2)]};
% % title(ax8,t8)
% %   ax8l=legend(ax8,'integer','prob','assigned');
%    ax8l=legend(ax8,t8);
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
    
    if exist('ax6','var')
        axt=ax6.copy;
    axt.Parent=f;
    subplot(3,3,8,axt) 
    end
    if exist('ax8','var')
        axt=ax8.copy;
    axt.Parent=f;
    subplot(3,3,9,axt) 
     ax8l=legend(axt,t8,'Location','northoutside');
%     legend(axt,'integer','prob','assigned','Location','northwestoutside');
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


pard.psfcheck.object=struct('String','PSF range','Style','checkbox');
pard.psfcheck.position=[3,1];

pard.PSFrange.object=struct('String','80 150','Style','edit');
pard.PSFrange.position=[3,2];
pard.PSFrange.Width=1;

pard.t3.object=struct('String','fit range histogram','Style','text');
pard.t3.position=[4,1];

pard.fitrange.object=struct('String','3 8','Style','edit');
pard.fitrange.position=[4,2];
pard.fitrange.Width=1;

pard.filecheck.object=struct('String','filenumbers','Style','checkbox');
pard.filecheck.position=[5,1];

pard.filenumbers.object=struct('String','1:100','Style','edit');
pard.filenumbers.position=[5,2];
pard.filenumbers.Width=1;

pard.plugininfo.type='ROI_Analyze';
  

pard.copy2page.object=struct('String','Copy to own page','Style','checkbox');
pard.copy2page.position=[6,1];
pard.copy2page.Width=2;

pard.plugininfo.type='ROI_Analyze';

end