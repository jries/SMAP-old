function stat=make_statistics2(locs,p,ploton)
if nargin<3
    ploton=true;
end

if p.filter
    modetxt={'layer','layer','layer','layer','layer','layer','layer','layer','layer','layer','layer','layer','layer','layer'};
else
    modetxt={'ungroup','group'};
end

if isfield(locs{1},'znm')&&~isempty(locs{1}.znm)
    txt='znm';
   zexist=true;
else
    txt='PSFxnm';
        zexist=false;
end

% fn=locD.files.file.name;
if ploton
if p.overview
    figure(34);


    sf=0.75;
    ax1=subplot(3,2,1);ax1.Position(3)=ax1.Position(3)*sf;
   
    ax2=subplot(3,2,2);ax2.Position(3)=ax2.Position(3)*sf;
    ax3=subplot(3,2,3);ax3.Position(3)=ax3.Position(3)*sf;
    ax4=subplot(3,2,4);ax4.Position(3)=ax4.Position(3)*sf;
    ax5=subplot(3,2,5);ax5.Position(3)=ax5.Position(3)*sf;
    ax6=subplot(3,2,6);ax6.Position(3)=ax6.Position(3)*sf;
else
    ax1=initaxis(p.resultstabgroup,'Photons');
    ax1.Position(3)=0.55;
    ax2=initaxis(p.resultstabgroup,'locprec');
    ax2.Position(3)=0.55;
    ax3=initaxis(p.resultstabgroup,'lifetime');
    ax3.Position(3)=0.55;
    ax4=initaxis(p.resultstabgroup,'BG');
    ax4.Position(3)=0.55;
    ax5=initaxis(p.resultstabgroup,txt);
    ax5.Position(3)=0.55;
    ax6=initaxis(p.resultstabgroup,['error ' txt]);
    ax6.Position(3)=0.55;
    if zexist
        ax7=initaxis(p.resultstabgroup,['error ' txt ' vs ' txt]);
        ax7.Position(3)=0.55;
    end


end
else
    ax1=[];ax2=[];ax3=[];ax4=[];ax5=[];ax6=[];ax7=[];
end
datrange=1:length(locs);


% look at frames
for k=datrange
    frames=locs{k}.frame;
    mf=max(frames);
    [hfr,n]=hist(frames,10);
    hfrc=hfr;
    hfrc(1:2)=[]; %ignore beginning
%     [~,ind]=max(diff(1./(hfr+.1*mean(hfr))));
%     ind=find(hfrc<(max(hfrc)-min(hfrc))/3+min(hfrc),1,'first');
    ind=find(hfrc<(max(hfrc))/3,1,'first');
    if isempty(ind)
        ind=length(n);
    else
    ind=ind+2;
    end
    r2=round([n(ind)-mf*.2 n(ind)+mf*.2]);
    frames2=frames(frames>=r2(1)&frames<=r2(2));
    [hfr2,n2]=hist(frames2,20);
    indco2=find(hfr2>(max(hfr2)-min(hfr2))/2+min(hfr2),1,'last');
    if isempty(indco2)
        indco2=length(n2);
    end
    falloffframe=n2(indco2);
    stat.frames.falloff(k)=falloffframe;
    [stat.frames.histogram(k).h,stat.frames.histogram(k).n]=hist(frames,100);
end

if ploton
    axf=initaxis(p.resultstabgroup,'frames');
    hold off
    for k=datrange
        plot(axf,stat.frames.histogram(k).n,stat.frames.histogram(k).h)
        hold on
        plot(axf,ones(2,1)*stat.frames.falloff(k),[0,max(stat.frames.histogram(k).h)])
    end
end

%photon stats
phot=getFieldAsVector(locs,'phot');
if isempty(phot{1})
    errdlg('no localizations in selected region')
    error('no localizations in selected region')
end
if p.checkphot
    for k=datrange
        phot{k}(phot{k}<p.photrange(1))=[];
        if length(p.photrange)>1
             phot{k}(phot{k}>p.photrange(2))=[];
        end
    end
    pr=p.photrange;
else
    pr=0.99;
end
hphot=plothist(phot,pr,[],0,ax1);
sphot={'Photons'};
phot1=1000;
phot2=3000;
for k=datrange
    sphot{end+1}='';
    sphot{end+1}=[num2str(k) '.' modetxt{k} ];
    Nloc(k)=length(phot{k});
    meanphot(k)=mean(phot{k});
    N1(k)=sum(phot{k}>phot1);
    N2(k)=sum(phot{k}>phot2);
    
    sphot{end+1}=['N'  ' = ' num2str(Nloc(k)/1000,'%5.0f') 'k'];
    sphot{end+1}=['<P'  '> = ' num2str(meanphot(k),'%5.0f')];
    sphot{end+1}=['r'  ' = ' num2str(N1(k)/N2(k),'%5.2f')];
    dat(k)=fitexpphot(hphot{k},[],ploton);
    sphot{end+1}=(['\mu'  ' = ' num2str(dat(k).mu,'%5.0f')]);   
end
stat.photons.Nloc=Nloc;
stat.photons.meanphot=meanphot;
stat.photons.mu=[dat(:).mu];

%locprec
locp=getFieldAsVector(locs,'locprecnm');
hlocp=plothist(locp,0.99,.25,0,ax2);
slp={'locprec_x'};
for k=datrange
    slp{end+1}='';
    slp{end+1}=[num2str(k) '.' modetxt{k} ];
    loch=locp{k};
    loch(loch<=0)=[];
    loch(loch>10000)=[];
    px = mylognfit(loch);
    [~,ind]=max(hlocp{k}.h);
    smx=hlocp{k}.n(ind);
    stat.locprec.max(k)=smx;
    slp{end+1}=['max: ' num2str(smx,3)];
    hf=mylognpdf(hlocp{k}.n,px(1),px(2))*sum(hlocp{k}.h)*(hlocp{k}.n(2)-hlocp{k}.n(1));
    if ploton
    plot(hlocp{k}.n,hf/max(hlocp{k}.h),'k:')
    end
    slp{end+1}=['median: ' num2str(median(locp{k}),3)];
    stat.locprec.median(k)=median(locp{k});
    indrise=find(hlocp{k}.h>1/2,1,'first');
    risingedge=hlocp{k}.n(indrise);
    stat.locprec.rising(k)=risingedge;
    slp{end+1}=['rising: ' num2str(risingedge,3)];
end

%lifetime
lifetime=getFieldAsVector(locs,'numberInGroup');
hlifet=plothist(lifetime,0.999,1,0,ax3);
slt={'lifetime'};
for k=datrange
    slt{end+1}='';
    slt{end+1}=[num2str(k) '.' modetxt{k} ];
    dat(k)=fitexpphot(hlifet{k},2,ploton);
    slt{end+1}=(['\mu'  ' = ' num2str(dat(k).mu,3)]);
    stat.lifetime.mu(k)=dat(k).mu;
end

%background
bg=getFieldAsVector(locs,'bg');
hbg=plothist(bg,0.95,1,0,ax4);
slb={'Background'};
for k=datrange
    slb{end+1}='';
    slb{end+1}=[num2str(k) '.' modetxt{k} ];
    mbg=mean(bg{k});
    slb{end+1}=['BG: ' num2str(mbg,'%5.0f')];
    stat.background.mean(k)=mbg;
end

%z/sigma
if zexist
    v=getFieldAsVector(locs,'znm');
else
    v=getFieldAsVector(locs,'PSFxnm');
end
hz=plothist(v,.99,[],0,ax5);
sls={txt};
for k=1:length(datrange)
    sls{end+1}='';
    sls{end+1}=[num2str(datrange(k)) '.' modetxt{datrange(k)} ];
    [~,ind]=max(hz{datrange(k)}.h);
    mx=hz{datrange(k)}.n(ind);    
    sls{end+1}=['max: ' num2str(mx,3)];
    stat.(txt).max(k)=mx;
end


if ploton
fontsize=14;
pos=[.7,0.025,.3,.95];
uicontrol('Parent',ax1.Parent,'style','text','String',sphot,'Units','normalized','Position',pos,'FontSize',fontsize,'HorizontalAlignment','left')
uicontrol('Parent',ax2.Parent,'style','text','String',slp,'Units','normalized','Position',pos,'FontSize',fontsize,'HorizontalAlignment','left')
uicontrol('Parent',ax3.Parent,'style','text','String',slt,'Units','normalized','Position',pos,'FontSize',fontsize,'HorizontalAlignment','left')
uicontrol('Parent',ax4.Parent,'style','text','String',slb,'Units','normalized','Position',pos,'FontSize',fontsize,'HorizontalAlignment','left')
uicontrol('Parent',ax5.Parent,'style','text','String',sls,'Units','normalized','Position',pos,'FontSize',fontsize,'HorizontalAlignment','left')
end

if zexist
    v=getFieldAsVector(locs,'locprecznm');
    hz=plothist(v,.99,[],0,ax6);
    slp={'locprecznm'};
    for k=datrange
        slp{end+1}='';
        slp{end+1}=[num2str(k) '.' modetxt{k} ];
        [~,ind]=max(hz{k}.h);
        mx=hz{k}.n(ind);    
        slp{end+1}=['max: ' num2str(mx,3)];
        stat.locprecznm.max(k)=mx;
    end    
    znm=getFieldAsVector(locs,'znm');
    rz=[-800 800];
    rsz=[0 100];
    him=myhist2(znm{1},v{1},10,1,rz,rsz);    
    if ploton
        axes(ax7)
        imagesc(rz,rsz,him')
        axis xy
        xlabel('znm');ylabel('locprec z');
    end
    if ploton
        uicontrol('Parent',ax6.Parent,'style','text','String',slp,'Units','normalized','Position',pos,'FontSize',fontsize,'HorizontalAlignment','left')
    end   
end

if ploton && ~p.overview
   ax1.Parent.Parent.SelectedTab=ax1.Parent;
end

function [v,datrange]=getvals(locD,field,p,indin)
if p.filter %use filtered values
    for layer=1:length(p.sr_layerson)
        if p.sr_layerson(layer)
            if p.useroi                
                v{layer}=locD.getloc(field,'layer',layer,'position','roi','within',indin).(field);
            else
                v{layer}=locD.getloc(field,'layer',layer,'within',indin).(field);
            end
            
        else
            v{layer}=0;
        end
    end
    datrange=find(p.sr_layerson);
else %use all values, plot for unconnected and connected
    if p.useroi
        position='roi';
    else
        position='all';
    end
       struc=locD.getloc(field,'position',position,'grouping','ungrouped','within',indin);
       v{1}=struc.(field);
       struc=locD.getloc(field,'position',position,'grouping','grouped','within',indin);
       v{2}=struc.(field);
    datrange=1:2;
end

function his=plothist(v,quantile,dphot,hmin,ax)

for k=1:length(v)
    if length(quantile)==1
    qq=myquantile(v{k},[1-quantile,quantile]);
    else
    qq=quantile;
    end
    q(k)=qq(2);q0(k)=qq(1);

    l(k)=length(v{k});
end

qmax=(max(q));
qmin=min(q0);
if qmax==qmin
    qmax=qmin+1;
end
lmax=max(l);
qfac=log10(lmax)-1;

if nargin==2||isempty(dphot)
dphot=(10^ceil(log10(qmax/qfac)))/100;
end
if nargin<4
    hmin=qmin;
end
nphot=hmin:dphot:qmax;
slegend={};
if ~isempty(ax)
    axes(ax)
    hold off
end
for k=1:length(v)
    
    if q(k)>0
%         sum(v{k})
        h=hist(v{k},nphot);
        [mmax,mi]=max(h(2:end-1)); 
        his{k}.h=h(2:end-1)/mmax;
        his{k}.n=nphot(2:end-1);
        if ~isempty(ax)
        plot(nphot(2:end-1),h(2:end-1)/mmax)
        hold on
        end
        slegend{end+1}=num2str(k);
    end
end
if ~isempty(ax)
legend(slegend,'Location','northeast')
end

function dat=fitexpphot(hin,fitstart,ploton)
h=double(hin.h);
xout=double(hin.n);
if length(h)>1
    [mmax,mi]=max(h(1:end-1)); 
    halft=find(h(mi:end)<mmax/2,1,'first')+mi;
if isempty(halft)
    halft=ceil(length(h)/2);
end
if nargin<2||isempty(fitstart)
    fitstart=ceil(mi*1.2);
end
fitr=fitstart:min(halft*5,length(h));

options=optimset('lsqcurvefit');
options.Display='off';
pf=lsqcurvefit(@expforfit,[1,xout(halft)],xout(fitr),h(fitr)/mmax,[],[],options);
if ploton
    plot(xout(fitr),expforfit(pf,xout(fitr)),'k--')
end
dat.mu=pf(2);
else
    dat.mu=0;  
end

function out=expforfit(p,x)
        out=p(1)*exp(-x/p(2));

