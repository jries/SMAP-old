function out=mywaveletfilter(in,level)
persistent wf1 wf2
if isempty(wf1)
    wf1=load('near_sym_a.mat');
    wf2=load('qshift_a.mat');
end
if nargin<2
    level=3;
end

[L,S]=dtwavexfm2_L(in,level,wf1,wf2);
out=dtwaveifm2_L(L,level,wf1,wf2,S);

if numel(out)~=numel(in)
    s=size(in);
    out=out(1:s(1),1:s(2));
end