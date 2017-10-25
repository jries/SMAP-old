global se
frametime=se.locData.files.file(1).info.timediff/1000; %in s
sites=se.sites;
dfall=[];
vall=[];
lenall=[];
timeall=[];
for k=1:length(sites)
    try
        polyh=sites(k).evaluation.BALM_fibril_growth.poly;
        dpol=polyh(2:end,:)-polyh(1:end-1,:);
    catch err
        continue
    end
    l=size(dpol,1);
    dfall(end+1:end+l)=dpol(:,2)*frametime;
     lenall(end+1:end+l)=dpol(:,1)/1000;
    vall(end+1:end+l)=dpol(:,1)./dpol(:,2)/frametime/1000; %in um/s
    timeall(end+1:end+l)=polyh(1:end-1,2)*frametime;
end

maxwaitelongation=0.1; %in um
stop=lenall<maxwaitelongation;


figure(99)
subplot(2,2,1);
histogram(dfall(stop),0:100:max(dfall(stop)));
title('wait time')
xlabel('wait time s')
subplot(2,2,2);
histogram(dfall(~stop),0:100:max(dfall(~stop)));
title('growth time')
xlabel('growth time s')


subplot(2,2,3);
histogram(lenall(~stop),0:.010:max(lenall(~stop)));
title('elongationstep')
xlabel('elongationstep um')

subplot(2,2,4);
histogram(vall(~stop),0:1e-3:max(vall(~stop)));
title('elongation speed')
xlabel('elongation speed um/s')

figure(100);
plot(timeall(stop),dfall(stop),'+')



