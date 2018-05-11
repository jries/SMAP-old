% challengescript
% parameters
dx=0; %corrections from bead fit.
dy=0;
dz=0;
photonfactor=1;
%density
densitysize_xy=40;
densitysize_z=80;
densitycutoff=3;
zmin=-600;zmax=600;
bgmin=15; bgmax=50; %background filter

ld=g.locData;

% group_dX=[15 30 50] 
locprec_cutoff=[20 40 30]; %for N1, N2, N3

PSF_cutoff_min=75; %only for Gaussian, so not relevant?
PSF_cutoff_max=200;
LLrel_cutoff=-2;

filenumber =1;

[~,filename]=fileparts(ld.files.file(filenumber).name);
if contains(filename,'.HD')
    density=1;
else
    density=0;
end
ind=strfind(filename,'.N');
if ~isempty(ind)
    brightness=str2num(filename(ind+2:ind+2));
else
    brightness=0;
end

group_dX=2*median(ld.loc.locprecnm);
if density==0
    group_dT=1;
else
    group_dT=0;
end
% determine HD, LD, bright, dark
%regroup with smaller dx, dt (depending on brightness, density)
ld.regroup(group_dX,group_dT)

%filter
indgood=true(size(ld.loc.xnm));
%  combined PSF grouping filter. PSF not relevant?
indgroupfilter=((...
    ld.loc.locprecnm<locprec_cutoff(brightness) ...
    & ld.loc.PSFxnm>PSF_cutoff_min & ld.loc.PSFxnm<PSF_cutoff_max ...
    & ld.loc.bg>bgmin & ld.loc.bg<bgmax  ...
    )| ld.loc.numberInGroup>1);

disp(['groupfilter: ' num2str(1-sum(indgroupfilter)/length(indgroupfilter))])
indgood=indgood & indgroupfilter;


%  LL filtering
if isfield(ld.loc,'LLrel') %do this filtr before regrouping?
    indgood=indgood & ld.loc.LLrel>LLrel_cutoff;
end

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
disp('done')
% write  dx, dy, dz from bead fit into challenge compare plugin

