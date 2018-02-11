% use rotation and shift from TOM to average localization data.
% for now: in ROImanager. Later maybe coordinate based

posfile='/Users/jonas/Documents/Data/Vilma/positions.txt';
em_motl_file='/Users/jonas/Documents/Data/Vilma/averaging/JR-Nup133_8fold_motl_11.em';

positions=csvread(posfile);
em=emread(em_motl_file);

iparticle=4;
idX=11;idY=12;idZ=13;
idX=14;idY=15;idZ=16;
iPhi=17;iPsi=18;iTheta=19;

pixrec=7; %(nm)

sites=g.locData.SE.sites;
%%
locnew=g.locData.loc;
            x0=nanmedian(locnew.xnm);
            y0=nanmedian(locnew.ynm);
newfile=g.locData.files.filenumberEnd+1;
used=false(size(locnew.xnm));
p.name='average';
p.addfile=true;


hold off
for siten=1:length(sites)
    g.status(num2str(siten));drawnow
siteh=sites(siten);
locsh1=g.locData.getloc({'xnm','ynm','znm'},'layer',1,'Position',siteh);
locsh2=g.locData.getloc({'xnm','ynm','znm'},'layer',2,'Position',siteh);

%make average and add as file
[locsite,indsite]=g.locData.getloc({'xnm','ynm','znm'},'layer',[1 2 ],'Position',siteh,'grouping','ungrouped');
used=used|indsite;

%transform to relative coordinates
clear c1 c2 ca;
c1(:,1)=locsh1.xnm-positions(siten,1);
c1(:,2)=locsh1.ynm-positions(siten,2);
c1(:,3)=locsh1.znm-positions(siten,3);

c2(:,1)=locsh2.xnm-positions(siten,1);
c2(:,2)=locsh2.ynm-positions(siten,2);
c2(:,3)=locsh2.znm-positions(siten,3);

ca(:,1)=locsite.xnm-positions(siten,1);
ca(:,2)=locsite.ynm-positions(siten,2);
ca(:,3)=locsite.znm-positions(siten,3);

% figure(88);
% scatter3(locsh1.xnmr,locsh1.ynmr,locsh1.znmr);
% hold on
% scatter3(locsh2.xnmr,locsh2.ynmr,locsh2.znmr)
% hold on
dc=[em(idX,siten) ,em(idY,siten),em(idZ,siten)]*pixrec*0;
angles=-[em(iPhi,siten),em(iPsi,siten),em(iTheta,siten)];
r1=applyT(c1,dc,angles);
% scatter3(r1x,r1y,r1z)

r2=applyT(c2,dc,angles);
% scatter3(r1(:,1),r1(:,2),r1(:,3))
% hold on
% scatter3(r2(:,1),r2(:,2),r2(:,3))

rall=applyT(ca,dc,angles);
locnew.xnm(indsite)=rall(:,1);locnew.ynm(indsite)=rall(:,2);locnew.znm(indsite)=rall(:,3);
locnew.filenumber(indsite)=newfile;
end
hold off


% from averageSites:
p.addfile=true;
fn=fieldnames(locnew);
%       for k=1:length(fn)
%            locnew.(fn{k})=locnew.(fn{k})(used);
%       end
     locc=g.locData.copy;
     locc.loc=locnew;
     locc.regroup;
     locc.filter;

    if p.addfile
            g.locData.addfile([p.name '_' num2str(newfile)]);
    initGuiAfterLoad(g);
    g.locData.SE.processors.preview.updateFilelist;
        locnew.xnm=locnew.xnm+x0;
        locnew.ynm=locnew.ynm+y0;
        for k=1:length(fn)
            g.locData.addloc(fn{k},locnew.(fn{k})(used))
        end
        g.locData.regroup;
        g.locData.filter;
    end

%%
function r=applyT(ri,dc, angles)
% ri(:,1)=ri(:,1)+dc(1);
% ri(:,2)=ri(:,2)+dc(2);
% ri(:,3)=ri(:,3)+dc(3);
r = tom_pointrotate(ri,angles(1),angles(2),angles(3));
r(:,1)=r(:,1)+dc(1);
r(:,2)=r(:,2)+dc(2);
r(:,3)=r(:,3)+dc(3);

end