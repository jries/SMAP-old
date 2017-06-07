%mammalian cells: temporal reconstruction
%parameters
p.cutoff=4; % for isosurface, at max(V(:))/cutoff
p.cutoffdc=3; % for isosurface, at max(V(:))/cutoff
p.sgauss=2;%smoothing of volume
p.sgauss2=.6;%smoothing of volume
p.winsize=200; %plotting
p.minx=0; %x: use 0 for central cut
p.pixelsize=5; %nm, for volume reconstruction
% p.limxy=[0 40]; % for plotting, in volume pixels
% p.limz=[0 40];
rot=1; % rotate sites to align bottom holes

timewindows=20; % number of sites per average: (number of sites / this number)
timepoints=100; % increment: (number of sites / this number)
framerate=25;

global se
sites=se.sites;


% rSphere=getFieldAsVector(sites,'evaluation','CME3DDSpherefit','map3D','rSphere');
mainFraction=getFieldAsVector(sites,'evaluation','CME3DDSpherefit','map3D','mainFraction');

[~,indsort]=sort(-mainFraction);

sitessort=sites(indsort);
% sitessort(1:10)=[];
rSphere=getFieldAsVector(sitessort,'evaluation','CME3DDSpherefit','map3D','rSphere');
mainFraction=getFieldAsVector(sitessort,'evaluation','CME3DDSpherefit','map3D','mainFraction');
topcoverage=getFieldAsVector(sitessort,'evaluation','CME3DDSpherefit','map3D','topcoverage');
thetac=getFieldAsVector(sitessort,'evaluation','CME3DDSpherefit','map3D','thetac');
phic=getFieldAsVector(sitessort,'evaluation','CME3DDSpherefit','map3D','phic');
thetacoverage=getFieldAsVector(sitessort,'evaluation','CME3DDSpherefit','map3D','bwfit','thetacoverage');

dN=round(length(indsort)/timewindows);
df=round(length(indsort)/timepoints);

range=1:df:length(indsort)-dN+df;
range(end)=min(range(end),length(sitessort)-dN);

rmean=zeros(1,length(indsort));
rcorr=zeros(1,length(indsort));
F(length(range))=struct('cdata',[],'colormap',[]);
for k=1:length(range)
    numloc=0;
    numloc2=0;
    
    indh=range(k):range(k)+dN;
    
    for s=indh
        numloc=numloc+length(sitessort(s).evaluation.CME3DDSpherefit.map3D.coordinates.r);  
        numloc2=numloc2+length(sitessort(s).evaluation.CME3DDSpherefit.map3D.coordinates2.x); 
    end
    
        
    thcov=mean(pi-thetacoverage(indh));
    Rhere=mean(rSphere(indh));
        rmean(indh)=Rhere;

    rcorr(indh)=rmean(indh)./rSphere(indh);   

    x=zeros(1,numloc);y=zeros(1,numloc);z=zeros(1,numloc);
    x2=zeros(1,numloc2);y2=zeros(1,numloc2);z2=zeros(1,numloc2);
    
    ind=1;ind2=1;
    for s=indh
%         if topcoverage(s)
%             tcorr=0;
%         else
%             tcorr=pi;
%         end
%         pcorr=0;
% %         tcorr=-thetac(s);pcorr=-phic(s);
        if rot
                xh=sitessort(s).evaluation.CME3DDSpherefit.map3D.coordinates.xc;
                x(ind:ind+length(xh)-1)=xh*rcorr(s);
                y(ind:ind+length(xh)-1)=sitessort(s).evaluation.CME3DDSpherefit.map3D.coordinates.yc*rcorr(s);
                z(ind:ind+length(xh)-1)=sitessort(s).evaluation.CME3DDSpherefit.map3D.coordinates.zc*rcorr(s);

                    xh2=sitessort(s).evaluation.CME3DDSpherefit.map3D.coordinates2.xc;
                    x2(ind2:ind2+length(xh2)-1)=xh2*rcorr(s);
                    y2(ind2:ind2+length(xh2)-1)=sitessort(s).evaluation.CME3DDSpherefit.map3D.coordinates2.yc*rcorr(s);
                    z2(ind2:ind2+length(xh2)-1)=sitessort(s).evaluation.CME3DDSpherefit.map3D.coordinates2.zc*rcorr(s);
        %         x(ind:ind+length(r)-1)=r;
        %         rch(ind:ind+length(r)-1)=r*rcorr(s);
        %         th(ind:ind+length(r)-1)=sitessort(s).evaluation.CME3DDSpherefit.map3D.coordinates.t+tcorr;
        %         ph(ind:ind+length(r)-1)=sitessort(s).evaluation.CME3DDSpherefit.map3D.coordinates.p+pcorr;

        else
                xh=sitessort(s).evaluation.CME3DDSpherefit.map3D.coordinates.x;
                x(ind:ind+length(xh)-1)=xh*rcorr(s);
                y(ind:ind+length(xh)-1)=sitessort(s).evaluation.CME3DDSpherefit.map3D.coordinates.y*rcorr(s);
                z(ind:ind+length(xh)-1)=sitessort(s).evaluation.CME3DDSpherefit.map3D.coordinates.z*rcorr(s);
                    xh2=sitessort(s).evaluation.CME3DDSpherefit.map3D.coordinates2.x;
                    x2(ind2:ind2+length(xh2)-1)=xh2*rcorr(s);
                    y2(ind2:ind2+length(xh2)-1)=sitessort(s).evaluation.CME3DDSpherefit.map3D.coordinates2.y*rcorr(s);
                    z2(ind2:ind2+length(xh2)-1)=sitessort(s).evaluation.CME3DDSpherefit.map3D.coordinates2.z*rcorr(s);        
        end        
        ind=ind+length(xh);
        ind2=ind2+length(xh2);
        
    end
    f=plotcoords(x,y,z,Rhere,thcov,p,x2,y2,z2);
    F(k)=getframe(f);
    

    
end
path=fileparts(se.files(1).name);
[f,path]=uiputfile([path filesep 'movie.mp4']);
if f
    v=VideoWriter([path f],'MPEG-4');
    v.FrameRate=framerate;
    open(v)
    for k=1:length(F)
        writeVideo(v,F(k));
    end
    close(v);
end
%f=figure(123);movie(f,F,2)



function f=plotcoords(x,y,z,R,thcov,p,x2,y2,z2)
% figure(88)
% scatter3(x,y,z+R,[],z+R)
% axis equal
% [t,p,r]=cart2sph(x,y,z);
% figure(89)
% polarplot(p,r,'.')
dz=R*cos(thcov);
dz=R;
xrange=[-p.minx,p.winsize];
yrange=[-p.winsize,p.winsize];
zrange=[-p.winsize,p.winsize]-p.winsize/2;

xr=xrange(1)+p.pixelsize:p.pixelsize:xrange(2)-p.pixelsize;
yr=yrange(1)+p.pixelsize:p.pixelsize:yrange(2)-p.pixelsize;
zr=zrange(1)+p.pixelsize:p.pixelsize:zrange(2)-p.pixelsize;
V=myhist3(x,y,-(z+dz),p.pixelsize,xrange,yrange,zrange);
Vs=smooth3(V,'gaussian',[1 1 1]*(round(2*p.sgauss/2)*2+1),p.sgauss);


V2=myhist3(x2,y2,-(z2+dz),p.pixelsize,xrange,yrange,zrange);
Vs2=smooth3(V2,'gaussian',[1 1 1]*(round(2*p.sgauss/2)*2+1),p.sgauss2);

[X,Y,Z]=meshgrid(yr,xr,zr);

co=max(Vs(:))/p.cutoff;
codc=max(Vs(:))/p.cutoffdc;

f1=figure(90);
clf
hiso=patch(isosurface(X,Y,Z,Vs,co),'EdgeColor','none','FaceColor',[1,.75,.65]);
isonormals(X,Y,Z,Vs,hiso);

hcap=patch(isocaps(X,Y,Z,Vs,co),'FaceColor','interp','EdgeColor','none');
axis equal
view(-35,30)

xlim(yrange)
ylim(xrange)
zlim(zrange)

lightangle(45,30);
lighting gouraud

f2=figure(91);
clf
hiso=patch(isosurface(X,Y,Z,Vs,codc),'EdgeColor','none','FaceColor',[1,.75,.65],'FaceAlpha',1);
isonormals(X,Y,Z,Vs,hiso);
hold on
if 0
indin=x2>0;
scatter3(y2(indin),x2(indin),-z2(indin)-dz,2,'k')
else
% codc=max(V2(:))/6;
codc=.5;
hiso2=patch(isosurface(X,Y,Z,Vs2,codc),'EdgeColor','none','FaceColor',[0,.75,1]);
isonormals(X,Y,Z,Vs2,hiso2);
hcap2=patch(isocaps(X,Y,Z,Vs2,co),'FaceColor','interp','EdgeColor','none');
end

axis equal
view(-35,30)

xlim(yrange)
ylim(xrange)
zlim(zrange)

lightangle(45,30);
lighting gouraud

f=f1;
% indg=r~=0;
% % [z,x]=pol2cart(p(indg),r(indg));
% [z,x,y]=sph2cart(t(indg),p(indg),r(indg));
% zm=z+mean(r);
% figure(87);
% plot(x,-zm,'.')
% % xlim([-200 200])
% % ylim([-200 200])
% % posp.x=-zm;posp.y=x;
% ps=5;
% edgesx=-200:ps:200;
% edgesz=-100:ps:300;
% srim=histcounts2(x,zm,edgesx,edgesz);
% [Y,X]=meshgrid(rangez(1)+pixelsz/2:pixelsz:rangez(2),ranger(1)+pixelsx/2:pixelsx:ranger(2));

% srim=histrender(posp,[-200 200], [-200 200], 5, 5);
% figure(89)
% imagesc(srim')
% axis('equal')

% % subplot(2,3,1)
% 
% imn=srim./(X'+pixelsx/4);
% % imn=srim';
% imn2=[imn(:,end) imn(:,end:-1:2) imn];

% figure(88);polarplot(th,rch,'.')
% rlim([0,max(rmean)*1.1])
drawnow;
% waitforbuttonpress
% pause(.2)

end