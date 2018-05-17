% challengescript
% parameters
%% DHNPC 10nm

% 44.453 45.166 9.744

dx=44.453; %corrections from bead fit.
dy=45.166;
dz=9.744;
photonfactor=0.662;

%% MT1.N1.LD optimized

%density
densitysize_xy=25;
densitysize_z=45;
densitycutoff=3;
zmin=-700;zmax=700;
bgmin=85; bgmax=105; %background filter
locprec_cutoff=30;
locprecz_cutoff=50;
phot_cutoff=1200;
LLrel_cutoff=-1.5;
group_dT=0;
border=10; %distance from min/max: if fit did converge to border



compare=true; %open compare java
% %% fit time
% fn=g.locData.files.file.name;
% l=load(fn);
% fitt=l.saveloc.fitparameters.processfittime;
% disp(['fit time without loading: ' num2str(fitt,3) ' s']);
%%

ld=g.locData;
maxiter=max(ld.loc.iterations);
% group_dX=[15 30 50] 


PSF_cutoff_min=75; %only for Gaussian, so not relevant?
PSF_cutoff_max=200;


filenumber =1;

[~,filename]=fileparts(ld.files.file(filenumber).name);


group_dX=min(50,2*median(ld.loc.locprecnm));

% determine HD, LD, bright, dark
%regroup with smaller dx, dt (depending on brightness, density)
ld.regroup(group_dX,group_dT)

%filter
indgood=true(size(ld.loc.xnm));
%  combined PSF grouping filter. PSF not relevant?
indgroupfilter=((...
    ld.loc.locprecnm<locprec_cutoff...
    & ld.loc.locprecznm<locprecz_cutoff...
    & ld.loc.PSFxnm>PSF_cutoff_min & ld.loc.PSFxnm<PSF_cutoff_max ...
    & ld.loc.bg>bgmin & ld.loc.bg<bgmax  ...
    & ld.loc.phot>phot_cutoff ...
    )| ld.loc.numberInGroup>1);

disp(['groupfilter: ' num2str(1-sum(indgroupfilter)/length(indgroupfilter))])
indgood=indgood & indgroupfilter;


%  LL filtering
if isfield(ld.loc,'LLrel') %do this filtr before regrouping?
    indgood=indgood & ld.loc.LLrel>LLrel_cutoff;
end

%iterations
indgood=indgood&ld.loc.iterations<maxiter;

%stripe artifacts: can be reduced if integer numbers are removed
indint=round(ld.loc.xnm)-ld.loc.xnm == 0 | ...
    round(ld.loc.ynm)-ld.loc.ynm == 0 | ...
    round(ld.loc.znm)-ld.loc.znm == 0;
indgood=indgood & ~indint;

%border filtering
indborder=ld.loc.xnm>min(ld.loc.xnm)+border & ...
    ld.loc.xnm<max(ld.loc.xnm)-border & ...
    ld.loc.ynm>min(ld.loc.ynm)+border & ...
    ld.loc.ynm<max(ld.loc.ynm)-border & ...
    ld.loc.znm>min(ld.loc.znm)+border & ...
    ld.loc.znm<max(ld.loc.znm)-border;
indgood=indgood&indborder;


% only current filenumber:
indgood =indgood &ld.loc.filenumber==filenumber;

%idea: modified cluster density: dx,y,z depends on median local locprecnm
%  density calculation, removal of single localisations
fdcal=figure(233);
dcal=plugin('Analyze','cluster','density_calculator',fdcal,g.P);
dcal.attachLocData(ld);
dcal.makeGui;
p=dcal.getGuiParameters;
p.countwhat.Value=1;
p.countingsize_xy=densitysize_xy;
p.countingsize_z=densitysize_z;
dcal.setGuiParameters(p);
dcal.useind=indgood;
dcal.processgo;
indcluster=ld.loc.clusterdensity>densitycutoff;

disp(['density filter: ' num2str(1-sum(indcluster)/length(indcluster))])
indgood=indgood & indcluster;


%  edges: z min, z max: remove or leave?
if isfield(ld.loc,'znm')
indgood = indgood & (ld.loc.znm>zmin & ld.loc.znm < zmax);
end


g.locData.setloc('challengefiltered',single(indgood));
g.locData.regroup(group_dX,group_dT);


% write x,y,z from grouped to ungrouped
copygroupfields={'xnm','ynm','znm'};
ldc=ld.copy; %dont overwrite in SMaP
[gi,sorti]=sort(ldc.loc.groupindex);
[gig,sortg]=sort(ldc.grouploc.groupindex);
inds=1;indg=1;
for k=1:gi(end)
    while gi(inds)<k
        inds=inds+1;
    end
    while gig(indg)<k
        indg=indg+1;
    end
    inds2=inds;
    while gi(inds2)==k && inds2<length(gi)
        inds2=inds2+1;
    end
    ind2g=indg;
    while gig(ind2g)==k && inds2<length(gig)
        ind2g=ind2g+1;
    end    
    indc=inds:inds2-1;
    gi(indc);
    gig(indg);
    if ~isempty(indc)
        if any(gi(indc)~=gig(indg))
            display('inconsitency group index')
        end
        for l=1:length(copygroupfields)
            if isfield(ldc.loc,copygroupfields{l})
                ldc.loc.(copygroupfields{l})(sorti(indc))=ldc.grouploc.(copygroupfields{l})(sortg(indg));
            end
        end
    end
end

if compare
fchallege=figure(234);
cs=plugin('Analyze','other','CompareToGroundTruthChallenge',fchallege,g.P);
cs.attachLocData(ldc);
cs.makeGui;
p=cs.getGuiParameters;
p.onlyfiltered=0;
p.offsetxyz=[dx dy dz];
p.photonfactor=photonfactor;
p.shiftframe=0;
cs.setGuiParameters(p);
cs.processgo;
end
disp('done')

% write  dx, dy, dz from bead fit into challenge compare plugin

%% create best J vs RMS plot
GTfile='/Volumes/t2ries/projects/SMLMChallenge2018/T1_MT0.N1.LD/activations.csv';
% GTfile='/Volumes/t2ries/projects/SMLMChallenge2018/T2_MT0.N1.LD/activations.csv';

dat = csvread(GTfile,1,0);
excessnoise=2;
phot=dat(:,6)*.9/excessnoise;
bg=90; 

PSF0=100; a=100; 

%PSF(z)
lambda=600;
w0=2*PSF0;
z=dat(:,5);
zR=pi*w0^2/lambda;
wz=w0*sqrt(1+(z/zR).^2);
PSF=wz/2;

locprecnm=sqrt((PSF.^2+a^2/12)./phot.*(16/9+8*pi*(PSF.^2+a^2/12)*bg./phot/a^2));
lps=sort(locprecnm);
norm=(1:length(lps))';
lpsj=sqrt(cumsum(lps.^2)./norm);
figure(88);
subplot(1,2,2)
hold off
plot(norm/max(norm),lpsj,'.')
hold on
plot([0 1],[1 1]*lpsj(1))
ax=gca;
ax.YLim(1)=0;
ax.YLim(2)=quantile(lpsj,0.99);
xlabel('Recall')
ylabel('RMS (nm)')
title('best possible RMS vs recall')

subplot(1,2,1);
histogram(locprecnm)
xlabel('localization precision nm')
xlim([0 quantile(locprecnm,0.98)])
title('localization precision Mortenson')
