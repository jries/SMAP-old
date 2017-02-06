function [imout,shiftedstack]=registerPSF3D(imin,p)
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
avim=mean(smallim,4);
for k=1:numbeads
    shift(k,:)=get3Dcorrshift(avim,smallim(:,:,:,k));
end
xn=1:size(imin,1);yn=1:size(imin,2);zn=1:size(imin,3);
[Xq,Yq,Zq]=meshgrid(yn,xn,zn);
%later: include 2x upscaling here
meanim=zeros(size(Xq));
% shift=-shift;
for k=1:numbeads
    shiftedstack(:,:,:,k)=interp3(imin(:,:,:,k),Xq-shift(1),Yq-shift(2),Zq-shift(3),'cubic',0);
    meanim=shiftedstack(:,:,:,k)+meanim;
end
imout=meanim/numbeads;
% average
%shift of everything to average
%shift volume


end