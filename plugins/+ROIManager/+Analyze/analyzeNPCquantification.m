
global se

fields={'evaluation','NPCLabelingQuantify'};
numfoundint=getFieldAsVector(se.sites,fields{:},'numfoundint');
numfoundrat=getFieldAsVector(se.sites,fields{:},'numfoundrat');
numbercornerassined=getFieldAsVector(se.sites,fields{:},'numbercornerassined');

nb=0:8;
    figure(89);
    
    subplot(3,3,1)
    hold off
    hi=hist(numfoundint,nb);
    bar(nb,hi)
    hold on
    pf=fithist(hi,nb);
    title(['int: ' num2str(pf,2)])
    
    axis tight
    subplot(3,3,2)
    hold off
    hr=hist(numfoundrat,nb);
    bar(nb,hr)<
    hold on
    pf=fithist(hr,nb);
    title(['rat: ' num2str(pf,2)])
    
    axis tight
    subplot(3,3,3)
    hold off
    ha=hist(numbercornerassined,nb);
    bar(nb,ha)
    hold on
    pf=fithist(ha,nb);
    title(['assign: ' num2str(pf,2)])
   axis tight
   
    subplot(3,3,4)
    hold off
    hall=hi+hr+ha;
    bar(nb,hall)
    hold on
    pf=fithist(hall,nb);
    title(['all: ' num2str(pf,2)])
    axis tight
    
    
if isfield(se.sites(1).evaluation.NPCLabelingQuantify,'numcornersfiltered') %not from simulation
    numcorners=getFieldAsVector(se.sites,fields{:},'numcornersfiltered');
    numcornersunf=getFieldAsVector(se.sites,fields{:},'numcorners');
    
    
    subplot(3,3,5)
    hold off
    hnc=hist(numcornersunf,nb);
    bar(nb,hnc)
    hold on
    pf=fithist(hnc,nb);
%     
%     hh=clusterfromlabeling(nb,8,2,.4)*sum(hnc);
%     plot(nb,hh)
    title(['true: ' num2str(pf,2)])
    axis tight
    subplot(3,3,[7,8,9])
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

function pf=fithist(hi,n)

shi=sum(hi);
% hi=hi/shi;

corners=8;rings=2;
% x=(0:corners)';
x=n(3:end)';
% clusterfromlabeling(x,corners,rings,.5)
ft=fittype('a*clusterfromlabeling(x,corners,rings,p)','problem',{'corners','rings'});
f=fit(x,hi(3:end)',ft,'problem',{corners, rings},'Lower',[0 0.01],'Upper',[inf .99],'Start',[shi .3]);
plot(0:8,f(0:8),'-g')
plot(x,f(x),'-*r')
pf=f.p;
end