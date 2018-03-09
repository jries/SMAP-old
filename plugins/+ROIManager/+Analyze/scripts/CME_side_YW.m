global se
sites=se.sites;
x1q13all = [];
x2q13all = [];
z1q13all = [];
z2q13all = [];
zmeddall = [];
for k=1:length(sites)
    try
        x1q13=sites(k).evaluation.CME2CSide_yule2.x1q13;
        x2q13=sites(k).evaluation.CME2CSide_yule2.x2q13;
        z1q13=sites(k).evaluation.CME2CSide_yule2.z1q13;
        z2q13=sites(k).evaluation.CME2CSide_yule2.z2q13;
        zmedd=sites(k).evaluation.CME2CSide_yule2.zmd;
    catch err
        continue
    end
    x1q13all = [x1q13all; x1q13];
    x2q13all = [x2q13all; x2q13];
    z1q13all = [z1q13all; z1q13];
    z2q13all = [z2q13all; z2q13];
    zmeddall = [zmeddall; zmedd];
end

figure(500)
subplot(2,3,1);
histogram(stepWidthAll, 'Binwidth', 1);
title('Rate')
xlabel('nm/frame')
subplot(2,3,2);
histogram(stallTimeAll, 'Binwidth', 1);
title('Stall time')
xlabel('Time')


figure(5000)
plot(1:2160', zmeddall)