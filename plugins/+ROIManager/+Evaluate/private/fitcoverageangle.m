function pos=fitcoverageangle(bwim,startp)
dangle=pi/64;
tr=-pi/2:dangle:pi/2;
pr=-pi:dangle:pi;
[Theta,Phi]=meshgrid(tr,pr);

lsqnonlin(@coverageerr,startph,[],[],[],Theta,Phi,bwim)

im=coverage(-pi/4,pi/4,pi/16,Theta,Phi);
figure(80);imagesc(tr,pr,im);
end

function im=coverage(ax,az,thetamin,Theta,Phi)
[T1,P1]=rotateSphericalCoordinates(Theta,Phi,1,ax);
[T2,P2]=rotateSphericalCoordinates(T1,P1,3,az);
im=T2>thetamin-pi/2;
end

function err=coverageerr(fitp,Theta,Phi,bwim)
imtest=coverage(fitp(1),fitp(2),fitp(3),Theta,Phi);
err=imtest(:)-bwim(:);
end