function plotTimeProfiles4Pi
global pathh
[f,pathh]=uigetfile([pathh '*.dcimg']);
if ~f
    return
end

il=imageloaderAll([pathh f]);
images=il.getmanyimages([],'mat');
images=images-myquantilefast(images(:),0.05,1000000);
numf=size(images,3);
fg=figure(127);
fg.Position(3)=1700; fg.Position(1)=10;
ax=axes(fg);
ax.Position=[0.05, 0.15,0.9,.8];
h.ax=ax;
h.findmax=uicontrol('Style','checkbox','String','center max','Units','normalized','Position',[0.8,0.05,0.05,0.05],'Value',1);
h.slider=uicontrol('Style','slider','Units','normalized','Position',[0.05,0.05,0.7,0.05],'Callback',{@slider_callback,images,h},'Max',numf,'Value',1,'Min',1);
h.roisize=uicontrol('Style','edit','String','5','Units','normalized','Position',[0.85,0.05,0.05,0.05]);
h.plot=uicontrol('Style','pushbutton','String','plot','Units','normalized','Position',[0.9,0.05,0.05,0.05],'Callback',{@plot_callback,images,h});
slider_callback(h.slider,0,images,h)
end

function slider_callback(a,b,images,h)
fr=ceil(a.Value);
imagesc(h.ax,images(:,:,fr))
axis(h.ax,'equal')
end

function plot_callback(a,b,images,h)
sigma=2;
pp=round(h.ax.CurrentPoint(1,1:2));
fr=ceil(h.slider.Value);
if h.findmax.Value
    hs=fspecial('gaussian',5,sigma);
    region=filter2(hs,images(pp(2)-5:pp(2)+5,pp(1)-5:pp(1)+5,fr));
    [~,indm]=max(region(:));
    [x,y]=ind2sub(size(region),indm);
    pp(2)=pp(2)-6+x;
    pp(1)=pp(1)-6+y;
end

f=figure(128);
rois=round(str2double(h.roisize.String));
dr=round(rois/2);
r=-dr:-dr+rois;
int=squeeze(sum(sum(images(pp(2)+r,pp(1)+r,:),1),2));

plot(int)
hold on
plot(0,0,'.')
% ax=gca;
% ax.YLim(1)=0;
end