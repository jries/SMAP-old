global se
ac=0*se.sites(1).evaluation.NPCgeomtryQuantify.ac;
z=se.sites(2).evaluation.NPCgeomtryQuantify.z;
zpos=0;
dzpos=55;
for k=1:length(se.sites)
    if (se.sites(k).pos(3))<zpos+dzpos && (se.sites(k).pos(3))>zpos-dzpos
    ac=ac+se.sites(k).evaluation.NPCgeomtryQuantify.ac;
    end
end
figure(88);
l=ceil(length(ac)/2);
zh=z(end-l+2:end);
ach=ac(2:l);
hold off
plot(zh,(ach)/max(abs(ach)));
hold on
dach=diff(ach);
ddach=diff(dach);
dz=z(2)-z(1);
plot(zh(1:end-1)+dz/2,(dach)/max(abs(dach)));
% plot(zh(1:end-2)+dz,(ddach)/max(abs(ddach)));
% ac=getFieldAsVector(se.sites,{'evaluation','NPCgeomtryQuantify','ac'})