function fh=cmeRingplotresults(p,results)
maxdrfac=p.maxdrro;

N=results.N;
ro=results.ro;
rc=results.rc;
dr=results.dr;
sigma=results.sigma;
% rdensity=results.rdensity;
filenumber=results.filenumber;
filenumberrange=results.filenumberrange;

drrel=dr./ro;
% drrelmax1=drrel;drrelmax1(drrelmax1>1)=1;
drrelmax=drrel;drrelmax(drrelmax>maxdrfac)=maxdrfac;
drri=dr;drri(drri>ro)=ro(drri>ro);


fh=figure(90);
subplot(3,4,1)
h=myhistfit(N,30,'normal');
[~,ind]=max(h(2).YData);
maxX=h(2).XData(ind); 
title(strcat('max ',num2str(maxX)));
% title(strcat(num2str(pd.mu),'±',num2str(pd.sigma)));
xlabel('N')
subplot(3,4,5)
plot(N,ro,'+')
xlabel('N')
ylabel('ro, dr')
hold on
plot(N,drri,'o')
hold off
legend('outer ring radius','ring thickness')
subplot(3,4,6)
plot(N,drrelmax,'x')
xlabel('N')
ylabel('dr/rout')

subplot(3,4,7)
plot(ro,drrelmax,'+')
xlabel('rout')
ylabel('dr/rout')

subplot(3,4,10)
% plot(mean(cat(3,rdensity{:}),3))
rdist=results.sumrdensityn(2)-results.sumrdensityn(1);
rnh=(0:length(results.sumrdensity)-1)*rdist;
plot(rnh,results.sumrdensity/results.numsites)
xlabel('radius (nm)')
ylabel('density')
title('radial density')

subplot(3,4,2)
hori=histogram(ro);
binedges=hori.BinEdges;
[h,pd]=myhistfit(ro);
xlabel('outer radius')
title(strcat(num2str(pd.mu),'±',num2str(pd.sigma)));

subplot(3,4,3)
% histogram(drrelmax);
% [~,pd]=myhistfit(drrelmax);
h=histfit(drrelmax,20,'kernel');
% title(strcat(num2str(pd.mu),'±',num2str(pd.sigma)));
[~,ind]=max(h(2).YData);
maxX=h(2).XData(ind); 
title(strcat('max ',num2str(maxX)));
xlabel('dr/rout');

subplot(3,4,9)
plot(filenumber,N,'+')
hold on
plot(filenumberrange,results.Nfilemean,'ok')
plot(filenumberrange,results.Nfilemedian,'*k')
hold off
xlabel('filenumber');
ylabel('N');
legend('N','Nmean','Nmedian')
subplot(3,4,4)
% histogram(rc,hori.BinEdges)
% [h,pd]=myhistfit(ro,hori.BinEdges);
[~,pd]=myhistfit(rc);
title(strcat(num2str(pd.mu),'±',num2str(pd.sigma)));
xlabel('radius circfit')

subplot(3,4,8)
% histogram(sigma);
[h,pd]=myhistfit(sigma);
title(strcat(num2str(pd.mu),'±',num2str(pd.sigma)));
xlabel('sigma imfit')

axim=subplot(3,4,12);

% xs=rc;ys=hori.BinEdges;
hdt=datacursormode(fh);
set(hdt,'UpdateFcn',{@CMElabel,results.sitenames,results.images,axim})

subplot(3,4,11);
imagesc(results.sumimage/results.numsites)
title('average image')
end