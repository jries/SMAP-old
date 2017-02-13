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
subplot(1,2,1)
nl=0:5:100;
h1=histogram(lenc,nl);
h1.FaceColor=[1 1 1]/2;

hold on
h2=histogram(len,nl);
h2.FaceColor=[.5 .5 1];
hold off
xlabel('distance (nm)')
legend('projected distance','raw distance')
subplot(1,2,2)
nang=-90:10:90;
histogram(anglec,nang)
xlabel('angle relative to cell axis (deg)')
title(['medP=' num2str(median(lenc),3) ', medr=' num2str(median(len),3) ', avP=' num2str(mean(lenc),3) ', avr=' num2str(mean(len),3)])


figure(99);histogram(lenc,20)
xlabel('projected distance (nm)')