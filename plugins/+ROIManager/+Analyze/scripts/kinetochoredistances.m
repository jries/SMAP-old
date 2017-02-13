%script to evaluate kinetochore distributions
sites=g.locData.SE.sites;
cells=zeros(length(sites),1);
angle=cells;
cellangle=cells;
len=cells;
for k=1:length(sites)
    si=sites(k);
    cells(k)=si.info.cell;
   
    %angle(k)=si.annotation.line1.angle;
    cellangle(k)=si.annotation.line2.angle;
    %len(k)=si.annotation.line1.length;
    [lenh,ah]=pos2len(si.annotation.line1.pos);
    angle(k)=ah/pi*180;
    len(k)=lenh*1000;
end


c=unique(cells);

for k=1:length(c)
    thisc=find(cells==c(k));
    ca=cellangle(thisc);
    ind=find(ca~=0,1,'first');
    if isempty(ind)
        ind=1;
    end
    cah=ca(ind);
    cellangle(thisc)=cah;
end

angler=(angle-cellangle);
angler(angler<-180)=angler(angler<-180)+360;
angler(angler>180)=angler(angler>180)-360;

anglec=angler;
anglec(anglec<-90)=anglec(anglec<-90)+180;
anglec(anglec>90)=anglec(anglec>90)-180;
lenc=abs(cos(anglec/180*pi)).*len;

figure(39)
subplot(2,2,1)
hold off
nl=0:5:100;
h1=histogram(lenc,nl);
h1.FaceColor=[1 1 1]/2;
hold on
r1=fitdist(h1.Data(h1.Data>0),'rician');
plot(h1.BinEdges+h1.BinWidth*0.5,r1.pdf(h1.BinEdges+h1.BinWidth*0.5)*sum(h1.BinCounts)*h1.BinWidth,'r')
xlabel('distance (nm)')
% legend('projected distance','raw distance')
xlim([0 max(max(len),max(lenc))+h1.BinWidth])
title(['Projected: d=' num2str(r1.s,'%1.1f') ' , sig=' num2str(r1.sigma,'%1.1f') ])


subplot(2,2,2)
hold off
h2=histogram(len,nl);
h2.FaceColor=[.5 .5 1];
hold on
r2=fitdist(h2.Data(h2.Data>0),'rician');
plot(h2.BinEdges+h2.BinWidth*0.5,r2.pdf(h2.BinEdges+h2.BinWidth*0.5)*sum(h2.BinCounts)*h2.BinWidth,'b')
hold off
xlabel('distance (nm)')
% legend('projected distance','raw distance')
xlim([0 max(max(len),max(lenc))+h1.BinWidth])

title(['raw: d=' num2str(r2.s,'%1.1f') ' , sig=' num2str(r2.sigma,'%1.1f')])
subplot(2,2,3)
nang=-90:10:90;
histogram(anglec,nang)
xlabel('angle relative to cell axis (deg)')
title(['medP=' num2str(median(lenc),3) ', medr=' num2str(median(len),3) ', avP=' num2str(mean(lenc),3) ', avr=' num2str(mean(len),3)])

subplot(2,2,4)
dy=len.*sin(angle);dx=len.*cos(angle);
plot(dx,dy,'x',mean(dx),mean(dy),'*')
title('displacement vectors')
% 
% 
% figure(99);histogram(lenc,20)
% xlabel('projected distance (nm)')