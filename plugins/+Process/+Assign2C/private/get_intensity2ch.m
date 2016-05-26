function loco=get_intensity2ch(loc,p)

%parameters

if isnumeric(p.split_slope)
    slope=p.split_slope;
else
    slope=str2num(p.split_slope);
end
if length(slope)==1;
    slope(2)=slope(1);
end


if isnumeric(p.split_offset)
    offset=p.split_offset;
else
    offset=str2double(p.split_offset);
end

if isnumeric(p.split_edge)
    se=p.split_edge;
else
se=str2num(p.split_edge);
end

edge1=se(1);
if length(se)>1
    edge2=se(2);
else
    edge2=edge1;
end

if isnumeric(p.split_intmin)
    imin=p.split_intmin;
else
imin=str2num(p.split_intmin);
end
int1min=imin(1);
if length(imin)>1
    int2min=imin(2);
else
    int2min=int1min;
end

loco.channel=0*loc.channel+3;
int1=loc.(p.assignfield1.selection);
int2=loc.(p.assignfield2.selection);
% int1=loc.intA1;
% int2=loc.intB1;

m1=myquantilefast(int1(:),[0.01,0.98],1e5);m2=myquantilefast(int2(:),[0.01,0.98],1e5);
map=max(m1(2), m2(2));mip=min(m1(1),m2(1));

ps=ceil((map-mip)/250);
indgood=int1~=0&int2~=0;
img=myhist2(int1(indgood),int2(indgood),ps,ps,[mip map],[mip map]);
% img=myhist2(int1,int2,ps,ps,[mip map],[mip map]);

ax1=initaxis(p.resultstabgroup,'i1 vs i2');
imagesc([mip map],[mip map],img)
hold on
plotboundary
hold off

ax2=initaxis(p.resultstabgroup,'log i1 vs i2');
limg=log(img);
limg(isinf(limg))=-1;
imagesc([mip map],[mip map],limg)
hold on
plotboundary
hold off

i1co=slope(2)*int2+edge1+offset(1);
c1=int1>i1co&int1>int1min;
loco.channel(c1)=1;


i2co=1./slope(1)*((int1)+edge2-offset(1));
c2=int2>i2co&int2>int2min;
loco.channel(c2)=2;

ax3=initaxis(p.resultstabgroup,'log split');
imgc2=myhist2(int1(c2),int2(c2),ps,ps,[mip map],[mip map]);
sout=size(imgc2);
outrgb=zeros(sout(1),sout(2),3);
outrgb(:,:,1)=imgc2;
imgc1=myhist2(int1(c1),int2(c1),ps,ps,[mip map],[mip map]);
outrgb(:,:,2)=imgc1;

c3=~(c1|c2);
imgc3=myhist2(int1(c3),int2(c3),ps,ps,[mip map],[mip map]);
outrgb(:,:,3)=imgc3;
outrgb=log(outrgb)+1;
outrgb(isinf(outrgb))=0;
imagesc([mip map],[mip map],outrgb/myquantilefast(outrgb(:),.9995,1e5))
hold on
plotboundary
hold off

function plotboundary
plot([mip map],slope(1)*[mip map]+offset(1),'w')
plot([mip map],slope(2)*[mip map]+offset(1),'w')
plot([mip map],slope(2)*[mip map]+edge1+offset(1),'m')
plot([mip map],slope(1)*[mip map]-edge2+offset(1),'m')
plot([mip 1/slope(2)*(int1min-edge1-offset(1))],(int1min)*[1 1],'c')
plot(int2min*[1 1],[mip slope(1)*int2min-edge2+offset(1)],'c')
end
end