function loco=get2Clocintensities(loc,transform,file,p)

 p.datapart.selection='all';
loct=apply_transform_locs(loc,transform,file,p);


[iA,iB,uiA,uiB]=matchlocsall(renamefields(loc),renamefields(loct),0,0,500);

loco.intA1=zeros(size(loc.xnm),'single');
loco.intB1=zeros(size(loc.xnm),'single');
loco.intA1(iA)=loc.phot(iA);
loco.intB1(iA)=loc.phot(iB);
loco.intA1(iB)=loc.phot(iA);
loco.intB1(iB)=loc.phot(iB);

function loco=renamefields(loci)
loco.x=loci.xnm;
loco.y=loci.ynm;
loco.frame=loci.frame;

