function out=myGPUwaveletfilter(in,level)
persistent H0 H1 H2 g1 g2
if isempty(H0)
    H0=3/8;
    H1=1/4;
    H2=1/16;
    g1=[H2,H1,H0,H1,H2];
    g2=[H2,0,H1,0,H0,0,H1,0,H2];
end

V1=conv2(conv2(in,g1','same'),g1,'same');
V2=conv2(conv2(V1,g2','same'),g2,'same');
out=im-V2;
end