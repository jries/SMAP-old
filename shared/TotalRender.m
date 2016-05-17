function [image,layers]=TotalRender(locData,pall,filterremove)
if nargin<3
    filterremove={};
end
imax=0;
if ~iscell(pall)
    pall={pall};
end
for k=1:length(pall)
    p=pall{k};
    if p.layercheck

        filterold=locData.getFilter(k);
        filternew=filterold;
        for f=1:length(filterremove)
            filternew=myrmfield(filternew,filterremove{f});
        end
        locData.setFilter(filternew,k);
        rawimage=renderSMAP(locData,p,k);
        locData.setFilter(filterold,k);
        layers(k).images.finalImages=drawerSMAP(rawimage,p);
        imax=max(imax,layers(k).images.finalImages.imax);

    end
end

image=displayerSMAP(layers,p);
image.imax=imax;
end