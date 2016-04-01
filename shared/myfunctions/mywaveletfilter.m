function out=mywaveletfilter(in,level,refine)
persistent wf1 wf2
if isempty(wf1)
    wf1=load('near_sym_a.mat');
    wf2=load('qshift_a.mat');
end
if nargin<2
    level=3;
end
if nargin<3
    refine=false;
end
% ino=in;
if refine
    level2=level*2;
    [L,S]=dtwavexfm2_L(in,level2,wf1,wf2);
    bg=dtwaveifm2_L(L,level2,wf1,wf2,S);
    in=in-bg;co=3*std(in(:)); in(in>co)=co;   
end
[L,S]=dtwavexfm2_L(in,level,wf1,wf2);
out=dtwaveifm2_L(L,level,wf1,wf2,S);
if refine
    out=out+bg;
end



if numel(out)~=numel(in)
    s=size(in);
    out=out(1:s(1),1:s(2));
end