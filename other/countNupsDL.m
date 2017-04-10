global path
[f,path]=uigetfile([path filesep '*.tif']);
%%
imgr=imageloaderAll([path f]);
img=imgr.getmanyimages(1:100,'mat');
%%

%%
imghr=imresize(img,1);

mimgr=max(imghr,[],3);

figure(39);imagesc(mimgr);
h=imfreehand;
mask=createMask(h);
mimg=double(mimgr).*double(mask);
figure(39)
imagesc(mimg);

maxima=maximumfindcall(mimg);
mint=maxima(:,3);
cutoff=600;
mintc=mint(mint>cutoff);
figure(33);
h=histogram(mintc,50);
cftool(h.BinEdges(1:end-1),h.Values)

maximac=maxima(mint>cutoff,:);
roih=3;
[X,Y]=meshgrid(-roih:roih);
intensity=[];stdx=[];stdy=[];as=[];
for k=size(maximac,1):-1:1
    x=maximac(k,2);y=maximac(k,1);
    if x>roih&& x < size(imgh,2)-roih &&y > roih &&y<size(imgh,1)-roih 
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

figure(33);
h=histogram(mintc(~badind),50);
cftool(h.BinEdges(1:end-1),h.Values)