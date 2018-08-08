function plotTimeProfiles4Pi
global pathh plotfig
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
h.findmax=uicontrol('Style','checkbox','String','center max','Units','normalized','Position',[0.825,0.05,0.05,0.05],'Value',1,'Parent',fg);
h.slider=uicontrol('Style','slider','Units','normalized','Position',[0.05,0.05,0.6,0.05],'Max',numf,'Value',1,'Min',1,'Parent',fg);
h.framenum=uicontrol('Style','edit','Units','normalized','Position',[0.65,0.05,0.025,0.05],'String','1','Parent',fg);


h.roisize=uicontrol('Style','edit','String','5','Units','normalized','Position',[0.875,0.05,0.025,0.05],'Parent',fg);


h.plotfft=uicontrol('Style','checkbox','Units','normalized','Position',[0.675,0.05,0.05,0.05],'String','FFT','Value',0,'Parent',fg);
h.fftrate=uicontrol('Style','edit','Units','normalized','Position',[0.7,0.05,0.025,0.05],'String','20','Parent',fg);
h.globalc=uicontrol('Style','checkbox','Units','normalized','Position',[0.725,0.05,0.05,0.05],'String','Glob contr','Parent',fg);
h.trendline=uicontrol('Style','checkbox','Units','normalized','Position',[0.775,0.05,0.05,0.05],'String','trend','Parent',fg);
h.plot=uicontrol('Style','pushbutton','String','plot','Units','normalized','Position',[0.9,0.05,0.05,0.05],'Callback',{@plot_callback,images,h},'Parent',fg);
h.slider.Callback={@slider_callback,images,h};
h.framenum.Callback={@slider_callback,images,h};
% plotfig=figure;
slider_callback(h.slider,0,images,h)
end

function slider_callback(a,b,images,h)
if strcmp(a.Style,'slider')
    fr=ceil(a.Value);
else
    fr=ceil(str2double(a.String));
end
implot=images(:,:,fr);
if h.globalc.Value
    implot(1,1)=max(images(:));
end
imagesc(h.ax,implot)
axis(h.ax,'equal')
h.framenum.String=num2str(fr);
h.slider.Value=fr;
end

function plot_callback(a,b,images,h)
global plotfig
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
if isempty(plotfig) || ~isvalid(plotfig)
    plotfig=figure;
end
f=figure(plotfig);
rois=round(str2double(h.roisize.String));
dr=round(rois/2);
r=-dr:-dr+rois;
intensity=squeeze(sum(sum(images(pp(2)+r,pp(1)+r,:),1),2));
 x=(1:length(intensity))';
plot(x,intensity)
hold on
plot(0,0,'.')
if h.trendline.Value
    fp=fit(x,intensity,'smoothingspline');
    plot(x,fp(x));
    
end


if h.plotfft.Value
    pwSpec=figure(129);
    Fs = str2double(h.fftrate.String);
    t = 0:1/Fs:(length(intensity)-1)*(1/Fs);

    N = length(intensity);

    xdft = fft(intensity);
    xdft = xdft(1:N/2+1);
    psdx = (1/(Fs*N)) * abs(xdft).^2;
    psdx(2:end-1) = 2*psdx(2:end-1);
    freq = 0:Fs/length(intensity):Fs/2;
    plot(freq,10*log10(psdx))
    grid on
    title('Periodogram Using FFT')
    xlabel('Frequency (Hz)')
    ylabel('Power/Frequency (dB/Hz)')
end
% ax=gca;
% ax.YLim(1)=0;
end