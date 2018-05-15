% challengescript
% parameters
%% MT1.N1.LD optimized
dx=0; %corrections from bead fit.
dy=0;
dz=0;
photonfactor=1;

%density
densitysize_xy=25;
densitysize_z=45;
densitycutoff=3;
zmin=-700;zmax=700;
bgmin=85; bgmax=105; %background filter
locprec_cutoff=20;
locprecz_cutoff=50;
phot_cutoff=200;
LLrel_cutoff=-2.5;
group_dT=0;
border=10; %distance from min/max: if fit did converge to border



compare=false; %open compare java
%%

ld=g.locData;
maxiter=max(ld.loc.iterations);
% group_dX=[15 30 50] 


PSF_cutoff_min=75; %only for Gaussian, so not relevant?
PSF_cutoff_max=200;


filenumber =1;

[~,filename]=fileparts(ld.files.file(filenumber).name);


group_dX=2*median(ld.loc.locprecnm);

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

copygroupfields={'xnm','ynm','znm'};





%  edges: z min, z max: remove or leave?
if isfield(ld.loc,'znm')
indgood = indgood & (ld.loc.znm>zmin & ld.loc.znm < zmax);
end


g.locData.setloc('challengefiltered',single(indgood));
g.locData.regroup(group_dX,group_dT);



% dcal=plugin('Analyze','cluster','density_calculator',fdcal,g.P);

% write x,y,z from grouped to ungrouped
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
cs.setGuiParameters(p);
cs.processgo;
end
disp('done')

% write  dx, dy, dz from bead fit into challenge compare plugin

