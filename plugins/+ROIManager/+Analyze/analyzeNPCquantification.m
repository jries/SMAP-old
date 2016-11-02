global se
fields={'evaluation','NPCLabelingQuantify'};
numcorners=getFieldAsVector(se.sites,fields{:},'numcornersfiltered');
numfoundint=getFieldAsVector(se.sites,fields{:},'numfoundint');
numfoundrat=getFieldAsVector(se.sites,fields{:},'numfoundrat');
numbercornerassined=getFieldAsVector(se.sites,fields{:},'numbercornerassined');

figure(89);
plot(numfoundint-numcorners)
hold on
plot(numfoundrat-numcorners)
plot(numbercornerassined-numcorners)
plot(numcorners/8);
hold off

legend('integer','prob','assigned','number of corners')