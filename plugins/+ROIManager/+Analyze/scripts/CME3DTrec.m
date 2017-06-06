%mammalian cells: temporal reconstruction
global se
sites=se.sites;

rot=true;
% rSphere=getFieldAsVector(sites,'evaluation','CME3DDSpherefit','map3D','rSphere');
mainFraction=getFieldAsVector(sites,'evaluation','CME3DDSpherefit','map3D','mainFraction');

[~,indsort]=sort(mainFraction);

sitessort=sites(indsort);
% sitessort(1:10)=[];
rSphere=getFieldAsVector(sitessort,'evaluation','CME3DDSpherefit','map3D','rSphere');
mainFraction=getFieldAsVector(sitessort,'evaluation','CME3DDSpherefit','map3D','mainFraction');
topcoverage=getFieldAsVector(sitessort,'evaluation','CME3DDSpherefit','map3D','topcoverage');
thetac=getFieldAsVector(sitessort,'evaluation','CME3DDSpherefit','map3D','thetac');
phic=getFieldAsVector(sitessort,'evaluation','CME3DDSpherefit','map3D','phic');

timewindows=10;
timepoints=50;
dN=round(length(indsort)/timewindows);
df=round(length(indsort)/timepoints);

range=1:df:length(indsort)-dN+df;
range(end)=min(range(end),length(sitessort)-dN);

rmean=zeros(1,length(indsort));
rcorr=zeros(1,length(indsort));
for k=1:length(range)
    numloc=0;
    
    indh=range(k):range(k)+dN;
    
    for s=indh
        numloc=numloc+length(sitessort(s).evaluation.CME3DDSpherefit.map3D.coordinates.r);  
    end
    
        

    Rhere=mean(rSphere(indh));
        rmean(indh)=Rhere;

    rcorr(indh)=rmean(indh)./rSphere(indh);
    
%     rh=zeros(1,numloc);rch=zeros(1,numloc);th=zeros(1,numloc);ph=zeros(1,numloc);
    x=zeros(1,numloc);y=zeros(1,numloc);z=zeros(1,numloc);
    ind=1;
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
%         x(ind:ind+length(r)-1)=r;
%         rch(ind:ind+length(r)-1)=r*rcorr(s);
%         th(ind:ind+length(r)-1)=sitessort(s).evaluation.CME3DDSpherefit.map3D.coordinates.t+tcorr;
%         ph(ind:ind+length(r)-1)=sitessort(s).evaluation.CME3DDSpherefit.map3D.coordinates.p+pcorr;

else
        xh=sitessort(s).evaluation.CME3DDSpherefit.map3D.coordinates.x;
        x(ind:ind+length(xh)-1)=xh*rcorr(s);
        y(ind:ind+length(xh)-1)=sitessort(s).evaluation.CME3DDSpherefit.map3D.coordinates.y*rcorr(s);
        z(ind:ind+length(xh)-1)=sitessort(s).evaluation.CME3DDSpherefit.map3D.coordinates.z*rcorr(s);
        
end        
        ind=ind+length(xh);
        
        
    end
    plotcoords(x,y,z,Rhere)
    

    
end


function plotcoords(x,y,z,R)
% figure(88)
% scatter3(x,y,z+R,[],z+R)
% axis equal
% [t,p,r]=cart2sph(x,y,z);
% figure(89)
% polarplot(p,r,'.')

V=myhist3(x,y,z+R,10,[-000,200],[-200,200],[-200,200]);
Vs=smooth3(V,'gaussian',[1 1 1]*3,1);

co=max(V(:))/8;
figure(90);
clf
hiso=patch(isosurface(Vs,co),'EdgeColor','none','FaceColor',[1,.75,.65]);
isonormals(Vs,hiso);

hcap=patch(isocaps(Vs,co),'FaceColor','interp','EdgeColor','none');
axis equal
view(35,30)

limvals=[0 50];
xlim(limvals)
ylim(limvals)
zlim(limvals)

lightangle(45,30);
lighting gouraud

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
pause(.2)

end