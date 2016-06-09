function [srim,nlocs,template]=constgaussrender(pos,rangex, rangey, pixelsx, pixelsy, lut,rangec,template)
[srim,nlocs,template]=histrender(pos,rangex, rangey, pixelsx, pixelsy, lut,rangec,template);
if isfield(pos,'s')
    s=myquantilefast(pos.s,0.5,1000)/pixelsx*ones(1,2);
    
elseif isfield(pos,'sx')
    s(2)=myquantilefast(pos.sx,0.5,1000)/pixelsx;
    s(1)=myquantilefast(pos.sy,0.5,1000)/pixelsy;
else
    s=[3 3];
end
fs=2*ceil(2.5*s)+1;
% s(2)=s(2)*10;
% h=fspecial('Gauss',ceil(5*s),s);
% s=size(srim);
% if length(s)==2
%     s(3)=1;
% end
% for k=1:s(3)
% srim(:,:,k)=filter2(h,srim(:,:,k));
% end

srim=imgaussfilt(srim,s,'FilterSize',fs);  
end