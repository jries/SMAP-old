function plotTimeProfiles4Pi
global pathh
[f,pathh]=uigetfile([pathh '*.dcimg']);
if ~f
    return
end

il=imageloaderAll([pathh f]);
images=il.getmanyimages([],'mat');
numf=size(images,3);
fg=figure(127);
fg.Position(3)=1700; fg.Position(1)=10;
ax=axes(fg);
ax.Position=[0.05, 0.15,0.9,.8];
h.ax=ax;
h.findmax=uicontrol('Style','checkbox','String','center max','Units','normalized','Position',[0.8,0.05,0.05,0.05]);
h.slider=uicontrol('Style','slider','Units','normalized','Position',[0.05,0.05,0.7,0.05],'Callback',{@slider_callback,images,h},'Max',numf);
h.roisize=uicontrol('Style','edit','String','3','Units','normalized','Position',[0.85,0.05,0.05,0.05]);
h.plot=uicontrol('Style','pushbutton','String','plot','Units','normalized','Position',[0.9,0.05,0.05,0.05],'Callback',{@plot_callback,images,h});

end

function slider_callback(a,b,images,h)
fr=ceil(a.Value);
imagesc(h.ax,images(:,:,fr))
axis(h.ax,'equal')
end

function plot_callback(a,b,images,h)

pp=round(h.ax.CurrentPoint(1,1:2));
fr=ceil(h.slider.Value);
if h.findmax.Value
    region=images(pp(2)-5:pp(2)+5,pp(1)-5:pp(1)+5,fr);
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
end