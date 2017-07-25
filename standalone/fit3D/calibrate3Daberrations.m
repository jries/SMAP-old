%  Copyright (c)2017 Ries Lab, European Molecular Biology Laboratory,
%  Heidelberg.
%  
%  This file is part of GPUmleFit_LM Fitter.
%  
%  GPUmleFit_LM Fitter is free software: you can redistribute it and/or modify
%  it under the terms of the GNU General Public License as published by
%  the Free Software Foundation, either version 3 of the License, or
%  (at your option) any later version.
%  
%  GPUmleFit_LM Fitter is distributed in the hope that it will be useful,
%  but WITHOUT ANY WARRANTY; without even the implied warranty of
%  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%  GNU General Public License for more details.
%  
%  You should have received a copy of the GNU General Public License
%  along with GPUmleFit_LM Fitter.  If not, see <http://www.gnu.org/licenses/>.
%  
%  
%  Additional permission under GNU GPL version 3 section 7
%  
%  If you modify this Program, or any covered work, by
%  linking or combining it with libraries required for interaction
%  with analysis programs such as Igor Pro or Matlab,
%  the licensors of this Program grant you additional permission
%  to convey the resulting work.
%%
function zcorr= calibrate3Daberrations(locs,pin)

%loc.filenumber
%loc.frame
%loc.phot
%loc.x
%loc.y
%loc.z
p.glassframe=[]; %automatic
p.setzero=true;
p.dz=10;
p.smoothz=1/1000;
p.cutoffrefine=150;
p.maxrange=800;
p=copyfields(p,pin);
if ~isfield(p,'smoothframe')
 p.smoothframe=2/10/p.dz;
end

f=figure('Name','Calibrate depht-induced aberrations');
p.tabgroup=uitabgroup(f);

%get beads from localizations
beads=segmentb_so(locs,p.dz);

%framerange
frange=[min(locs.frame) max(locs.frame)];

% get true positions f0 for beads
for k=length(beads):-1:1
    [beads(k).f0]=getf0Z_so(beads(k).loc,p.dz);
    beads(k).phot=beads(k).loc.phot(min(length(beads(k).loc.phot),max(1,round(beads(k).f0))));
end
f0all=([beads(:).f0]);
badind=f0all<frange(1)|isnan(f0all); 
beads(badind)=[];

% determine position of the glass: no beads below glass
if isempty(p.glassframe)
    p.axhere=axes(uitab(p.tabgroup,'Title','f0'));
    p.glassframe=getf0glass(beads,p);
end
p.f0glass=p.glassframe;% *ones(1,max([beads(:).filenumber]));

%calculate relevant other coordinates
axh=axes(uitab(p.tabgroup,'Title','z vs f0'));

hold off
for k=1:length(beads)
    beads(k).loc.zglass=(beads(k).loc.frame-p.f0glass)*p.dz;
    beads(k).loc.z0relative=-(beads(k).f0-beads(k).loc.frame)*p.dz;
    beads(k).loc.dzcorr=beads(k).loc.z0relative-beads(k).loc.z;
    beads(k).f0glass=beads(k).f0-p.f0glass;
    beads(k).loc.z0glass=beads(k).f0glass*p.dz+0*beads(k).loc.zglass;
    indplot=abs(beads(k).loc.z0relative)<p.maxrange;
    plot(axh,beads(k).loc.zglass(indplot),beads(k).loc.z(indplot),'.')
    hold on
end
xlabel('frame')
ylabel('fitted z position (nm)')
p.axhere=[];
phere=p;
phere.smoothing=[0.05 0.002];

% iteratively determine spline approximation and remove beads that
% are too far away
[ZcorrInterp]=getZinterp(beads,[],phere,'zglass');
phere=p;
phere.cutoffrefine=500;
[ZcorrInterp]=getZinterp(beads,ZcorrInterp,phere);
[err1,dzerr]=geterrors(beads,ZcorrInterp,p);
goodind=find(true(length(beads),1));
beads2=beads;
while  1% length(beads2)>length(beads)/2
    cutoff=3*nanmean(err1);
    badind=(err1>cutoff|isnan(err1));
    if sum(badind)==0
        break
    end
    goodind=goodind(~badind);
    beads2=beads(goodind);
    [ZcorrInterp]=getZinterp(beads2,ZcorrInterp,p);
    %calculate errors
    [err1,dzerrh]=geterrors(beads2,ZcorrInterp,p);    
    dzerr(goodind)=dzerrh;
end

%plot output
p.axhere=axes(uitab(p.tabgroup,'Title','Interpolation'));
[ZcorrInterp]=getZinterp(beads2,ZcorrInterp,p);

%correct beads for testing
ax1=axes(uitab(p.tabgroup,'Title','Validation'));
f=ax1.Parent;
ax2=axes(f,'Position',[0.5 0 1 1]);
subplot(1,2,1,ax1);
subplot(1,2,2,ax2);

minv=inf;
maxv=-inf;
for k=length(beads):-1:1
   dZ=ZcorrInterp.interp(beads(k).loc.zglass,beads(k).loc.z);
   beads(k).loc.zcorrected=beads(k).loc.z+dZ;
   if any(goodind==k)
       inr=abs(beads(k).loc.z0relative)<1000;
       if sum(inr)>0
           minv=min(min(minv,min(beads(k).loc.z(inr))),min(beads(k).loc.zcorrected(inr)));
           maxv=max(max(maxv,max(beads(k).loc.z(inr))),max(beads(k).loc.zcorrected(inr)));
       end
   else
       col='r.';
       plot(ax1,beads(k).loc.z0relative,beads(k).loc.z,col)
       plot(ax2,beads(k).loc.z0relative,beads(k).loc.zcorrected,col)
       hold(ax1,'on');
       hold(ax2,'on');
   end
end
for k=length(beads):-1:1
    if any(goodind==k)
       plot(ax1,beads(k).loc.z0relative,beads(k).loc.z,'k.')
       plot(ax2,beads(k).loc.z0relative,beads(k).loc.zcorrected,'k.')
       hold(ax1,'on');
       hold(ax2,'on');
    end
end

xlim(ax1,[-1000 1000])
ylim(ax1,[minv maxv]);
xlim(ax2,[-1000 1000])
ylim(ax2,[minv maxv]);

xlabel(ax1,'true z (nm)')
ylabel(ax1,'fitted z (nm)')

xlabel(ax2,'true z (nm)')
ylabel(ax2,'corrected z (nm)')

zcorr=ZcorrInterp.interp;    

end


function [Zint]=getZinterp(beads,Zintold,p,zaxis)
if nargin<4
    zaxis='zglass';
end
 %make big array for interpolation
 zglassall=[];z0relativeall=[];zfitall=[];idall=[];zplot=[];dzerrall=[];z0glassall=[];
for k=1:length(beads)
    zglassall=double(vertcat(zglassall,beads(k).loc.zglass));
    z0glassall=double(vertcat(z0glassall,beads(k).loc.z0glass));
    z0relativeall=double(vertcat(z0relativeall,beads(k).loc.z0relative));
    zfitall=double(vertcat(zfitall,beads(k).loc.z));
    idall=double(vertcat(idall,k*ones(length(beads(k).loc.zglass),1)));
    zplot=double(vertcat(zplot,beads(k).loc.dzcorr));          
end

if strcmp(zaxis,'z0glass')
    zzax=z0glassall;
else
    zzax=zglassall;
end

qzfit=myquantile(zfitall,[0.05,0.95]);
qzfit(1)=qzfit(1)+p.dz;qzfit(2)=qzfit(2)-p.dz;
inz=abs(z0relativeall)<p.maxrange;
inz=inz&(zfitall)<qzfit(2)&(zfitall)>qzfit(1);
inz=inz&abs(zplot)<p.maxrange;
if ~isempty(Zintold)  
    dz=Zintold.interp(zzax,zfitall)-zplot;
    inz=inz&abs(dz)<p.cutoffrefine;
    h=histcounts(idall(inz),(1:max(idall)+1))';
    minpoints=p.maxrange/p.dz;
    innump=h(idall)>minpoints;
    inz=inz&innump;
end
zfitx=(qzfit(1):p.dz:qzfit(2))';
zfitallh=vertcat(zfitall(inz),zfitx,zfitx,zfitx);
zploth=vertcat(zplot(inz),0*zfitx,0*zfitx,0*zfitx);
zzaxh=vertcat(zzax(inz),min(zzax)+0*zfitx,min(zfitx)+0*zfitx,mean([min(zfitx),min(zzax)])+0*zfitx);          

xrange=min(zzax(inz)):100:max(zzax(inz));

yrange=[qzfit(1):10: qzfit(2)];
[X,Y]=meshgrid(xrange,yrange);  
%interpolation
Z=RegularizeData3D(zzaxh,zfitallh,zploth,xrange,yrange,'smoothness',[p.smoothframe p.smoothz],'extend','always');
Zint.interp=griddedInterpolant(X',Y',Z');
Zint.zaxis=zaxis;

if ~isempty(p.axhere)
    Zplot=Zint.interp(X,Y);
    minzp=min(zfitall(inz));
    Zplot(Zplot<minzp)=minzp-10;
    error=abs(Zint.interp(zzax(inz),zfitall(inz))-zplot(inz));
    scatter3(p.axhere,zzax(inz),zfitall(inz),zplot(inz),5,(1-error/max(error))*min(Zplot(:)))
    xlabel(p.axhere,'objective position above glass (nm)');ylabel(p.axhere,'zfit (nm)'); zlabel(p.axhere,'correction (nm)');
    hold(p.axhere,'on')
    s=surf(p.axhere,X,Y,Zplot);
    s.FaceAlpha=0.8;
    s.EdgeColor='none';
    p.axhere.ZLim(1)=minzp;
end

end

function  f0glass=getf0glass(beads,p)
%determine the position of the glass as the robust minimum of bead
%positions
if isempty(p.axhere)
    f=figure;ax=gca;
else
    ax=p.axhere;
end

f0=[beads.f0];
dzh=50/p.dz;
induse=f0<dzh*60;
f0=f0(induse);    
range=min(f0):dzh:max(f0);
h=histogram(ax,f0,range);


[mh]=max(h.Values);
ind=find(h.Values>mh*.4,1,'first');
f0h=range(ind);
ind=(f0>f0h-2*dzh&f0<f0h+2*dzh);
f0glass=mean(f0(ind));
hold(ax, 'on');
xlabel('bead position (frame)')
ylabel('counts')
title('bead positions')

if isempty(p.axhere)
    close(f)
else
    plot(ax,f0glass,ones(size(f0glass)),'k*')
end
end

function [err1,dzerr]=geterrors(beads,Zint,p)
yrange=Zint.interp.GridVectors{2};
for k=1:length(beads)
    zh=double(beads(k).loc.z);
    zglass=beads(k).loc.zglass;
    z0f=beads(k).loc.dzcorr;
    inz=abs(zh<300) & abs(z0f)<300 & (zh)<yrange(end) & (zh)>yrange(1);
    dz=Zint.interp(zglass(inz),zh(inz))-z0f(inz);
    if p.setzero&&beads(k).f0glass*p.dz<2*p.dz %on glass
      factor=0.2;
    else
      factor=1;
    end        
    err1(k)=mean(dz.^2)*factor;
    err2(k)=mean(abs(dz))*factor;
    err3(k)=std(dz)*factor;
    dzerr{k}=Zint.interp(zglass,zh)-z0f;
end
end

