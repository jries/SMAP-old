function pos=fitcoverageangle(bwim,startp)
bwim=2*bwim-1;
sin=size(bwim);
% dangle(1)=2*pi/sin(1);
% dangle(2)=pi/sin(2);
% % dangle=pi/64;
% tr=-pi/2:dangle(2):pi/2;
% pr=-pi:dangle(1):pi;
tr=linspace(-pi/2,pi/2,sin(2));
pr=linspace(-pi,pi,sin(1));

[Theta,Phi]=meshgrid(tr,pr);

startph=[0.3,0.3,pi/2];
imstart=coverage(startph(1),startph(2),startph(3),Theta,Phi);
options=optimset('lsqnonlin');
options.Algorithm='levenberg-marquardt';
[fitp,resnorm]=lsqnonlin(@coverageerr,startph,[],[],options,Theta,Phi,bwim);
[fitp2,resnorm2]=lsqnonlin(@coverageerr,startph,[],[],options,Theta,Phi,-bwim);

if resnorm2<resnorm
    fitp=fitp2;
end
imfit=coverage(fitp(1),fitp(2),fitp(3),Theta,Phi);
% im=coverage(-pi/4,pi/4,pi/16,Theta,Phi);
imcombine=zeros(size(imstart,1),size(imstart,2),3);
imcombine(:,:,1)=bwim;
imcombine(:,:,2)=imfit;
imcombine(:,:,3)=imstart;
figure(80);imagesc(tr,pr,imcombine);
figure(81);imagesc(tr,pr,imfit);
end

function im=coverage(ax,az,thetamin,Theta,Phi)
[T1,P1]=rotateSphericalCoordinates(Theta,Phi,3,az);
[T2,P2]=rotateSphericalCoordinates(T1,P1,1,ax);
% im=2*((T2thetamin-pi/2)-1);
im=T2-(thetamin-pi/2);
co=pi/64;
im=im/co;
im(im<-1)=-1;im(im>1)=1;
end

function err=coverageerr(fitp,Theta,Phi,bwim)
% fitp
% fitp(3)=pi/4;
imtest=coverage(fitp(1),fitp(2),fitp(3),Theta,Phi);
err=imtest(:)-bwim(:);
if 1
imcombine=zeros(size(bwim,1),size(bwim,2),3);
imcombine(:,:,1)=bwim;
imcombine(:,:,2)=imtest;
% imcombine(:,:,3)=imstart;
figure(77);imagesc(imcombine);title(fitp)
drawnow
end
end