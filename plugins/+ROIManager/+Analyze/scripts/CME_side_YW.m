global se
sites=se.sites;
global sites
numSites = length(sites);
siteIdx = 1:30:numSites;
[siteIdx,~] = meshgrid(siteIdx, 1:30);
siteIdx = siteIdx(:)';
x1q13all = zeros(1, numSites);
x2q13all = [];
z1q13all = [];
z2q13all = [];
zmeddall = [];
dPeakAll = zeros(1, numSites);
N1all = zeros(1, numSites);
N2all = zeros(1, numSites);
for k=1:numSites
    try
        x1q13all(k)=sites(k).evaluation.CME2CSide_yule2.x1q13;
        x2q13=sites(k).evaluation.CME2CSide_yule2.x2q13;
        z1q13=sites(k).evaluation.CME2CSide_yule2.z1q13;
        z2q13=sites(k).evaluation.CME2CSide_yule2.z2q13;
        zmedd=sites(k).evaluation.CME2CSide_yule2.zmd;
        dPeak = max(sites(k).evaluation.CME2CSide_yule2.dPeak);
        if length(dPeak)==1
            dPeakAll(k) = dPeak;
        end
        N1all(k)=sites(k).evaluation.CME2CSide_yule2.N1;
        N2all(k)=sites(k).evaluation.CME2CSide_yule2.N2;
    catch err
        continue
    end
    x2q13all = [x2q13all; x2q13];
    z1q13all = [z1q13all; z1q13];
    z2q13all = [z2q13all; z2q13];
    zmeddall = [zmeddall; zmedd];
end



while 0
    figure(20400)
    opts = optimoptions('ga','PlotFcn',@gaplotbestf);
    [x, yval, fg]= ga(@findMinCor, 300, [repelem(1,300);repelem(-1,300)],[150;-150],[],[],repelem(0,300),repelem(1,300),[],1:300,opts);
    xy1G1 = getKernelMatrix(true, [500 500], [xnmrot1G1 ynmrot1G1]);
    xy2G1 = getKernelMatrix(true, [500 500], [xnmrot2G1 ynmrot2G1]);
    xy1G2 = getKernelMatrix(true, [500 500], [xnmrot1G2 ynmrot1G2]);
    xy2G2 = getKernelMatrix(true, [500 500], [xnmrot2G2 ynmrot2G2]);
    figure(2002)
    subplot(2,2,1);
    imagesc(xy1G1)
    subplot(2,2,2);
    imagesc(xy2G1)
    subplot(2,2,3);
    imagesc(xy1G2)
    subplot(2,2,4);
    imagesc(xy2G2)
    
    xy1D = xcorr2(xy1G1,xy1G2);
    max(xy1D(:))

    xy2D = xcorr2(xy2G1,xy2G2);
    max(xy2D(:))

end

xnmrot1Avg = [];
ynmrot1Avg = [];
xnmrot2Avg = [];
ynmrot2Avg = [];


for i = 1:numSites
    xnmrot1Avg = [xnmrot1Avg; sites(i).evaluation.CME2CSide_yule2.locs1.xnmrot];
    ynmrot1Avg = [ynmrot1Avg; sites(i).evaluation.CME2CSide_yule2.locs1.ynmrot];
    xnmrot2Avg = [xnmrot2Avg; sites(i).evaluation.CME2CSide_yule2.locs2.xnmrot];
    ynmrot2Avg = [ynmrot2Avg; sites(i).evaluation.CME2CSide_yule2.locs2.ynmrot];
end
avg = getKernelMatrix(true, [500 500], [xnmrot1Avg ynmrot1Avg; xnmrot2Avg ynmrot2Avg]);


Y = [];
for i = 1:numSites
    singleOne = getKernelMatrix(true, [500 500], [sites(i).evaluation.CME2CSide_yule2.locs1.xnmrot sites(i).evaluation.CME2CSide_yule2.locs1.ynmrot; sites(i).evaluation.CME2CSide_yule2.locs2.xnmrot sites(i).evaluation.CME2CSide_yule2.locs2.ynmrot]);
    D = xcorr2(avg, singleOne);
    Y(i) = log10(max(D(:)));
end
transY=(Y-std(Y))/mean(Y);

figure(500)
subplot(2,3,1);
scatter(siteIdx, dPeakAll);
title('Distance (density peaks)')

subplot(2,3,3);
scatter(siteIdx, x2q13all);
title('Interquartile range (x)')

subplot(2,3,4);
scatter(siteIdx, z2q13all);
title('Interquartile range (z)')

subplot(2,3,2);
scatter(siteIdx, -zmeddall);
title('Distance (median)')

subplot(2,3,5);
scatter(siteIdx, transY);
title('XCorr')


featureMatrix = [(transY'-1)*100 -zmeddall N1all' N2all'];
%featureMatrix = [transY' z2q13all -zmeddall N1all' N2all'];

tSNE = tsne(featureMatrix);
figure(501);
gscatter(tSNE(:,1),tSNE(:,2),siteIdx, '','',10);

wanderlustR=wanderlust(featureMatrix);
wdr = mean(wanderlustR.T);

figure(502)
scatter(siteIdx',wdr')

[~,IWdr] = sort(wdr');
pcorW = corr(siteIdx',siteIdx(IWdr)', 'Type', 'Spearman');
sum(siteIdx'==siteIdx(IWdr)')/numSites

[~,IZmedd] = sort(-zmeddall');
pcorZmedd = corr(siteIdx',siteIdx(IZmedd)', 'Type', 'Spearman');
sum(siteIdx'==siteIdx(IZmedd)')/numSites

TY = (transY'-1)*100;
[~,ITY] = sort(TY');
pcorTY = corr(siteIdx',siteIdx(ITY)', 'Type', 'Spearman');
sum(siteIdx'==siteIdx(ITY)')/numSites

function y = findMinCor(x)
    global sites
    xnmrot1G1 = [];
    ynmrot1G1 = [];
    xnmrot2G1 = [];
    ynmrot2G1 = [];
    xnmrot1G2 = [];
    ynmrot1G2 = [];
    xnmrot2G2 = [];
    ynmrot2G2 = [];

    idxG1=find(x);
    idxG2=find(~x);
    for i = idxG1
        xnmrot1G1 = [xnmrot1G1; sites(i).evaluation.CME2CSide_yule2.locs1.xnmrot];
        ynmrot1G1 = [ynmrot1G1; sites(i).evaluation.CME2CSide_yule2.locs1.ynmrot];
        xnmrot2G1 = [xnmrot2G1; sites(i).evaluation.CME2CSide_yule2.locs2.xnmrot];
        ynmrot2G1 = [ynmrot2G1; sites(i).evaluation.CME2CSide_yule2.locs2.ynmrot];
    end

    for i = idxG2
        xnmrot1G2 = [xnmrot1G2; sites(i).evaluation.CME2CSide_yule2.locs1.xnmrot];
        ynmrot1G2 = [ynmrot1G2; sites(i).evaluation.CME2CSide_yule2.locs1.ynmrot];
        xnmrot2G2 = [xnmrot2G2; sites(i).evaluation.CME2CSide_yule2.locs2.xnmrot];
        ynmrot2G2 = [ynmrot2G2; sites(i).evaluation.CME2CSide_yule2.locs2.ynmrot];
    end
    G1 = getKernelMatrix(true, [500 500], [xnmrot1G1 ynmrot1G1; xnmrot2G1 ynmrot2G1]);
    G2 = getKernelMatrix(true, [500 500], [xnmrot1G2 ynmrot1G2; xnmrot2G2 ynmrot2G2]);
    
    D = xcorr2(G1,G2);
    y=log10(max(D(:)));
end

function [z, bw] = getKernelMatrix(centralized,matrixSize, sVrot, varargin) 
    [meshx, meshy] = meshgrid(0:matrixSize(1), 0:matrixSize(2)); % make full grid
    if centralized
        shx = matrixSize(1)/2; shy = matrixSize(2)/2;
    else
        shx = 0; shy = 0;
    end
    q = [meshx(:)-shx, meshy(:)-shy]; % convert the grid to positions
    if length(varargin)==0
        [sVrotx,xy,bw]=ksdensity(sVrot, q); % get kernel density estimation
    else
        [sVrotx,xy,bw]=ksdensity(sVrot, q, varargin{1,1}, varargin{1,2}); % get kernel density estimation    
    end
    
    % converted into an image
    idx = sub2ind(matrixSize(end:-1:1)+1, xy(:,2)+1+shx, xy(:,1)+1+shy); 
    z = zeros(matrixSize(end:-1:1)+1);
    z(idx) = sVrotx;
end
