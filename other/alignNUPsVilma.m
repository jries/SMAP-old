% use rotation and shift from TOM to average localization data.
% for now: in ROImanager. Later maybe coordinate based

posfile='/Users/jonas/Documents/Data/Vilma/positions.txt';
em_motl_file='/Users/jonas/Documents/Data/Vilma/averaging/JR-Nup133_8fold_motl_11.em';
parIDfile='/Users/jonas/Documents/Data/Vilma/averaging/ParID.mat';

positions=csvread(posfile);
em=emread(em_motl_file);
l=load(parIDfile);
parID=l.ParID;
iparticle=4;
idX=11;idY=12;idZ=13;
% idX=14;idY=15;idZ=16;
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
p.addfile=false;

thetan=-pi:pi/64:pi;
actheta1=zeros(size(thetan));
actheta2=zeros(size(thetan));

zn=-100:2:100;
acz1=zeros(size(zn));
acz2=zeros(size(zn));

hold off
for k=1:length(parID)
    
    expr='_S(?<id>\d+)C';
    m=(regexp(parID{k,2},expr,'names'));
    ID=str2double(m.id);
    siten=g.locData.SE.indexFromID(sites,ID);
    
    
    g.status(num2str(siten));drawnow
siteh=sites(siten);
locsh1=g.locData.getloc({'xnm','ynm','znm'},'layer',1,'Position',siteh);
locsh2=g.locData.getloc({'xnm','ynm','znm'},'layer',2,'Position',siteh);

%make average and add as file


%transform to relative coordinates
clear c1 c2 ca;
c1(:,1)=locsh1.xnm-positions(siten,1);
c1(:,2)=locsh1.ynm-positions(siten,2);
c1(:,3)=locsh1.znm-positions(siten,3);

c2(:,1)=locsh2.xnm-positions(siten,1);
c2(:,2)=locsh2.ynm-positions(siten,2);
c2(:,3)=locsh2.znm-positions(siten,3);




% figure(88);
% scatter3(locsh1.xnmr,locsh1.ynmr,locsh1.znmr);
% hold on
% scatter3(locsh2.xnmr,locsh2.ynmr,locsh2.znmr)
% hold on
dc=-[em(idX,k) ,em(idY,k),em(idZ,k)]*pixrec;
angles=[em(iPhi,k),em(iPsi,k),em(iTheta,k)];
r1=applyT(c1,dc,angles);
% scatter3(r1x,r1y,r1z)

r2=applyT(c2,dc,angles);
% scatter3(r1(:,1),r1(:,2),r1(:,3))
% hold on
% scatter3(r2(:,1),r2(:,2),r2(:,3))
if p.addfile
    [locsite,indsite]=g.locData.getloc({'xnm','ynm','znm'},'layer',[1 2 ],'Position',siteh,'grouping','ungrouped');
    used=used|indsite;
    ca(:,1)=locsite.xnm-positions(siten,1);
    ca(:,2)=locsite.ynm-positions(siten,2);
    ca(:,3)=locsite.znm-positions(siten,3);
    rall=applyT(ca,dc,angles);
    % rall=ca;
    locnew.xnm(indsite)=rall(:,1);locnew.ynm(indsite)=rall(:,2);locnew.znm(indsite)=rall(:,3);
    locnew.filenumber(indsite)=newfile;
end


%analysis
theta=cart2pol(r1(:,1),r1(:,2));
histtheta=hist(theta,thetan);
histtheta2=hist(theta+pi,thetan+pi);
tac=myxcorr(histtheta,histtheta);
tac=tac+myxcorr(histtheta2,histtheta2);

actheta1=actheta1+tac;


histz=hist(r1(:,3),zn);
acz1=acz1+myxcorr(histz,histz);

theta=cart2pol(r2(:,1),r2(:,2));
histtheta=hist(theta,thetan);
histtheta2=hist(theta+pi,thetan+pi);
tac=myxcorr(histtheta,histtheta);
tac=tac+myxcorr(histtheta2,histtheta2);
actheta2=actheta2+tac;

histz=hist(r2(:,3),zn);
acz2=acz2+myxcorr(histz,histz);
end
hold off

normz=length(zn):-1:1;
normt=length(thetan):-1:1;
figure(109);
subplot(2,1,1);
plot(zn,acz1/acz1(2)./normz,zn,acz2/acz2(2)./normz);
subplot(2,1,2);
plot(thetan,actheta1/actheta1(2)./normt,thetan,actheta2/actheta2(2)./normt);

% from averageSites:
    if p.addfile
p.addfile=true;
fn=fieldnames(locnew);
%       for k=1:length(fn)
%            locnew.(fn{k})=locnew.(fn{k})(used);
%       end
     locc=g.locData.copy;
     locc.loc=locnew;
     locc.regroup;
     locc.filter;


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