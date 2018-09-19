sites=g.locData.SE.sites;
aca=0;
for k=1:length(sites)
    if isfield(sites(k).evaluation.NPCgeomtryQuantify,'actheta')
    ach=sites(k).evaluation.NPCgeomtryQuantify.actheta;
    aca=ach+aca;
    end
end

tn=sites(1).evaluation.NPCgeomtryQuantify.thetan;
figure(88);hold on;
norm=length(tn)-(1:length(tn));
plot(tn(2:end)*50,aca(2:end)./norm(2:end));