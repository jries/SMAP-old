function [imout,shiftedstack,shift]=registerPSF3D(imin,p)
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
smallim=zeros(length(p.xrange),length(p.yrange),length(p.framerange),size(imin,4));
for k=1:size(imin,4)
    smallim(:,:,:,k)=imin(p.xrange,p.yrange,p.framerange,k);
end
avim=mean(imin,4);

xn=1:size(imin,1);yn=1:size(imin,2);zn=1:size(imin,3);
[Xq,Yq,Zq]=meshgrid(yn,xn,zn);
meanim=zeros(size(Xq));
refim=avim(p.xrange,p.yrange,p.framerange);
for k=numbeads:-1:1
%     shift(k,:)=get3Dcorrshift(avim,smallim(:,:,:,k));
    shift(k,:)=get3Dcorrshift(refim,smallim(:,:,:,k));
    shiftedstack(:,:,:,k)=interp3(imin(:,:,:,k),Xq-shift(k,2),Yq-shift(k,1),Zq-shift(k,3),'cubic',0);
    meanim=shiftedstack(:,:,:,k)+meanim;
    refim=meanim(p.xrange,p.yrange,p.framerange);
end

%later: include 2x upscaling here

%  shift=-shift;
for k=numbeads:-1:1
    
    
end
imout=meanim/numbeads;

f=figure(99);
delete(f.Children)
% avim(2,2,2)=max(avim(:));
imageslicer(vertcat(avim,imout,shiftedstack(:,:,:,1),imin(:,:,:,1)),f);

% average
%shift of everything to average
%shift volumea


end