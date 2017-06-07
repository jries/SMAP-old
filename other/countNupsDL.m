global path
[f,path]=uigetfile([path filesep '*.tif']);
%%
imgr=imageloaderAll([path f]);
imga=double(imgr.getmanyimages(1:100,'mat'));
%
offset=10000;
cutoffmin=200;
img=imga-offset;
quantile(img(:),0.02)
%%
imghr=imresize(img,1);

mimgr=max(imghr,[],3);

figure(39);imagesc(mimgr);
h=imfreehand;
mask=createMask(h);
mimg=double(mimgr).*double(mask);
figure(39)
imagesc(mimg);
mmmm=mimg(mimg~=0);
background=quantile(mmmm(:),0.2);
mimg=mimg-background;
%%
numbins=100; 
maxima=maximumfindcall(mimg);
mint=maxima(:,3);
% cutoffmin=600;
mintc=mint(mint>cutoffmin);
figure(33);
hold off
h=histogram(mintc,numbins);
% cftool(h.BinEdges(1:end-1)+h.BinWidth/2,h.Values)
xfit=h.BinEdges(1:end-1)+h.BinWidth/2;
yfit=h.Values;
fitp=fit(xfit',yfit','gauss1','Robust','LAR');
hold on
plot(xfit,fitp(xfit))
title(['maximum at ' num2str(fitp.b1,4) ', bg ' num2str(background)])
% cftool(h.BinEdges(1:end-1),h.Values)

%%

maximac=maxima(mint>cutoffmin,:);
roih=3;
[X,Y]=meshgrid(-roih:roih);
intensity=[];stdx=[];stdy=[];as=[];
for k=size(maximac,1):-1:1
    x=maximac(k,2);y=maximac(k,1);
    if x>roih&& x < size(imghr,2)-roih &&y > roih &&y<size(imghr,1)-roih 
        imsmall=double(mimg(x-roih:x+roih,y-roih:y+roih));
        imsmall=imsmall-min(imsmall(:));
        [as(k),alpha(k)]=asymmetry(imsmall);
        mx=sum(sum(imsmall.*X))/sum(imsmall(:));
        my=sum(sum(imsmall.*Y))/sum(imsmall(:));
        stdx(k)= sum(sum(imsmall.*(X-mx).^2))/sum(imsmall(:));
        stdy(k)= sum(sum(imsmall.*(Y-my).^2))/sum(imsmall(:));
        imrs=imresize(imsmall,8);
        intensity(k)=max(imrs(:));
    end
end



figure(34);
subplot(2,2,1)
plot(as,intensity,'.');
xlabel('asymmetry');
ylabel('intensity');

subplot(2,2,2);
plot(stdx,intensity,'.');

xlabel('stdx');
ylabel('intensity');

subplot(2,2,3);
plot(stdy,intensity,'.');

xlabel('stdy');
ylabel('intensity');

subplot(2,2,4);
plot(stdx./stdy,intensity,'.');

xlabel('stdx/stdy');
ylabel('intensity');

ascutoff=0.15; stdcutoff=3.2;
badind=as>ascutoff|stdx>stdcutoff|stdy>stdcutoff;
badind=badind+mintc'<cutoffmin;

figure(35);
hold off
h=histogram(mintc(~badind),numbins);
% cftool(h.BinEdges(1:end-1),h.Values)

% cftool(h.BinEdges(1:end-1)+h.BinWidth/2,h.Values)
xfit=h.BinEdges(1:end-1)+h.BinWidth/2;
yfit=h.Values;
fitp=fit(xfit',yfit','gauss1','Robust','LAR');
hold on
plot(xfit,fitp(xfit))
title(['Filtered. Maximum at ' num2str(fitp.b1,5) ', bg ' num2str(background)])
