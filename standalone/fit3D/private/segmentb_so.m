function bead=segmentb_so(locs,dz)
% if isfield(p,'zrangeuse')
% minframes=(p.zrangeuse(2)-p.zrangeuse(1))/ p.dz/2;
% else
    minframes=800/dz;
% end
% minframes=(p.zrangeuse(2)-p.zrangeuse(1))/ p.dz/5;

% locDatacopy=obj.locData.copy;
% locDatacopy.regroup(150,minframes/2);

sortmatrix=horzcat(locs.filenumber,locs.frame,locs.x);
[~,indsort]=sortrows(sortmatrix,1:3);
list=connectsingle2c(double(locs.x(indsort)),double(locs.y(indsort)),double(locs.frame(indsort)),double(0.5),int32(0),int32(10000));
hh=histcounts(list,1:max(list)+1);
beadind=find(hh>=minframes);
% beadind=find(locDatacopy.grouploc.numberInGroup>minframes);
% xb=locDatacopy.grouploc.xnm(beadind);
% yb=locDatacopy.grouploc.ynm(beadind);
% filenumber=locDatacopy.grouploc.filenumber(beadind);
% locg=locDatacopy.getloc({'xnm','ynm','numberInGroup','filenumber'},'layer',1,'Position','roi','removeFilter','filenumber','grouping','grouped');
% beadind=find(locg.numberInGroup>minframes);
% xb=locg.xnm(beadind);
% yb=locg.ynm(beadind);
% filenumber=locg.filenumber(beadind);
winsize=2.5;
for k=length(beadind):-1:1
    indh=list==beadind(k);indh1=find(indh,1,'first');
    beadposx=median(locs.x(indh));
    beadposy=median(locs.y(indh));
    beadfile=locs.filenumber(indh1);
%     beadfile=locs.filenumber(indh1);
    inrange=(locs.x-beadposx).^2+(locs.y-beadposy).^2<winsize^2&locs.filenumber==beadfile;
    %     beadhere=(mywithin(locDatacopy.loc.xnm,[xb(k)-winsize/2 winsize],locDatacopy.loc.ynm,[yb(k)-winsize/2 winsize]));
%     beadfile=locDatacopy.loc.filenumber==filenumber(k);
%     beadhere=beadhere&beadfile;
%     bead(k).loc=locDatacopy.getloc({'xnm','ynm','PSFxnm','PSFynm','frame','phot','znm'},'filenumber',filenumber(k),'Position',[xb(k) yb(k) winsize]);
%     bead(k).loc=copystructReduce(locDatacopy.loc,beadhere,{'xnm','ynm','PSFxnm,','PSFynm','frame'});
    bead(k).filenumber=beadfile;
    bead(k).pos=[beadposx beadposy];
    bead(k).loc.z=double(locs.z(inrange));
    bead(k).loc.phot=double(locs.phot(inrange));
    bead(k).loc.frame=double(locs.frame(inrange)); 
end

%remove beads with close neighbours.
mindist=10;
goodind=true(length(beadind),1);
for k=1:length(beadind)
    for l=k+1:length(beadind)
        d2=sum((bead(k).pos-bead(l).pos).^2);
%         df=mean(bead(k).loc.frame)-mean(bead(l).loc.frame);
        if bead(k).filenumber==bead(l).filenumber && d2<mindist^2 %&& df<p.zrangeuse(2)-p.zrangeuse(1)
            goodind(k)=false;
            goodind(l)=false;
        end
    end
end
bead=bead(goodind);
end
