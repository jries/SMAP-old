global se
sites=se.sites;
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
    catch err
        continue
    end
    x2q13all = [x2q13all; x2q13];
    z1q13all = [z1q13all; z1q13];
    z2q13all = [z2q13all; z2q13];
    zmeddall = [zmeddall; zmedd];
end

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

featureMatrix = [x1q13all' x2q13all z2q13all -zmeddall];

tSNE = tsne(featureMatrix);
figure(501);
gscatter(tSNE(:,1),tSNE(:,2),siteIdx, '','',10);


