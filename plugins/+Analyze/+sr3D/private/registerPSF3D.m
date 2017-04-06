function [imout,shiftedstackn,shift,indgood]=registerPSF3D(imin,p,axs)
if nargin<3
    axs={};
end
% perform correlation on:
% p.xrange
% p.yrange
% p.framerange
% can do 2D if length(p.framerange)=1
%cutout small volumes
if ~isfield(p,'xrange')
    p.xrange=1:size(imin,1);
end
if ~isfield(p,'yrange')
    p.yrange=1:size(imin,2);
end
if ~isfield(p,'framerange')
    p.framerange=1:size(imin,3);
end

numbeads=size(imin,4);
if numbeads==1
    imout=imin;
    shiftedstack=imin;
    shift=[0 0 0 ];
    indgood=true;
    return
end
smallim=zeros(length(p.xrange),length(p.yrange),length(p.framerange),size(imin,4));
for k=1:size(imin,4)
    frh=round(p.framerange-p.zshiftf0(k));
    smallim(:,:,:,k)=imin(p.xrange,p.yrange,frh,k);
end
avim=nanmean(imin,4);

xn=1:size(imin,1);yn=1:size(imin,2);zn=1:size(imin,3);
[Xq,Yq,Zq]=meshgrid(yn,xn,zn);
% meanim=zeros(size(Xq));
meanim=[];
refim=avim(p.xrange,p.yrange,p.framerange);

simin=size(imin);
shiftedstack=zeros(simin(1),simin(2),simin(3),numbeads)+NaN;
if ~isempty(p.zshiftf0)
    zshiftf0=p.zshiftf0;
else
    zshiftf0=zeros(numbeads,1);
end
for k=1:numbeads
%     shift(k,:)=get3Dcorrshift(avim,smallim(:,:,:,k));
    if p.alignz
        [shift(k,:),cc(k)]=get3Dcorrshift(refim,smallim(:,:,:,k));
    else
        [shift(k,:),cc(k)]=get2Dcorrshift(refim,smallim(:,:,:,k));
    end
    
    shiftedh=interp3(imin(:,:,:,k),Xq-shift(k,2),Yq-shift(k,1),Zq-shift(k,3)-double(zshiftf0(k)),'cubic',0);
%     if ~isempty(meanim)
%         ratio=shiftedh(p.xrange,p.yrange,p.framerange)./meanim(p.xrange,p.yrange,p.framerange);
%         factor=nanmedian(ratio(:));
%     else
%         factor=nanmax(shiftedh(:));
%     end
%     if factor==0
%         factor=NaN;
%     end
    shiftedstack(:,:,:,k)=shiftedh;%/factor;
%     shiftedh(isnan(shiftedh))=0;
%     meanim=shiftedh+meanim;
    meanim=nanmean(shiftedstack(:,:,:,1:k),4);
    
    refim=meanim(p.xrange,p.yrange,p.framerange);
end

shiftedstackn=normalizstack(shiftedstack,p);


[indgood]=getoverlap(shiftedstackn,cc,shift,p,true(1,size(shiftedstackn,4)));
[indgood]=getoverlap(shiftedstackn,cc,shift,p,indgood);
[indgood,res,normglobal,co]=getoverlap(shiftedstackn,cc,shift,p,indgood);

shiftedstackn=shiftedstackn/normglobal;

imout=nanmean(shiftedstackn(:,:,:,indgood),4);
shiftedstackn(1,end,:,~indgood)=nanmax(shiftedstackn(:));
shiftedstackn(1,:,1,~indgood)=nanmax(shiftedstackn(:));

if length(axs)>0
hold(axs{1},'off')
plot(axs{1},(res(~indgood)),cc(~indgood),'rx');title(co)
hold(axs{1},'on')
plot(axs{1},(res(indgood)),cc(indgood),'g*');title(co)
xlabel(axs{1},'residulas')
ylabel(axs{1},'cross-correlation value')
end

if length(axs)>1
    imageslicer(vertcat(avim,imout),'Parent',axs{2}.Parent)
end
%later: include 2x upscaling here

%  shift=-shift;
% for k=numbeads:-1:1
%     
%     
% end
% imout=meanim/numbeads;

% f=figure(99);
% delete(f.Children)
% % avim(2,2,2)=max(avim(:));
% imageslicer(vertcat(avim,imout,shiftedstack(:,:,:,1),imin(:,:,:,1)),'Parent',f);
% imageslicer(shiftedstack);
% average
%shift of everything to average
%shift volumea


end

function [indgood,res,normamp,co]=getoverlap(shiftedstackn,cc,shift,p,indgood)
refimn=nanmean(shiftedstackn(p.xrange,p.yrange,p.framerange,indgood),4);
normamp=nanmax(refimn(:));
shiftedstackn=shiftedstackn/normamp;
refimn=refimn/normamp;
for k=size(shiftedstackn,4):-1:1
     sim=shiftedstackn(p.xrange(2:end-1),p.yrange(2:end-1),p.framerange,k);
     dv=(refimn(2:end-1,2:end-1,:)-sim).^2;
    res(k)=sqrt(nansum(dv(:)));
    sumall(k)=nansum(sim(:));
end
rescc=res./cc;
rescc(abs(shift(:,1))>3|abs(shift(:,2))>3)=NaN;
[a,b]=robustMean(rescc(cc>0));
if isnan(b)
    a=nanmean(rescc);b=nanstd(rescc);
end
co=a+3.*b;
indgood=rescc<co;
end

function out=normalizstack(in,p)
sin=size(in);
out=0*in+NaN;
midp=round((length(p.xrange)+1)/2);
xr=p.xrange(midp-3:midp+3);yr=p.yrange(midp-3:midp+3);
for k=1:sin(4)
    imh=in(xr,yr,p.framerange,k);
    nm=nanmean(imh(:));
    if nm>0
    out(:,:,:,k)=in(:,:,:,k)/nm;
    end
end
% meanim=nanmean(out,4);

% out=0*in+NaN;
% for k=1:sin(4)
%     imh=in(:,:,:,k);
%     ratio=imh(p.xrange,p.yrange,p.framerange)./meanim(p.xrange,p.yrange,p.framerange);
%     factor=nanmedian(ratio(:));
%     if factor>0
%     out(:,:,:,k)=imh/factor;
%     end
% end

% iii=squeeze(out(round(sin(1)/2),round(sin(1)/2),:,:));
% figure(88);plot(iii)
end