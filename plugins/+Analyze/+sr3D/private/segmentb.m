function bead=segmentb(obj,p)

minframes=(p.zrangeuse(2)-p.zrangeuse(1))/ p.dz/2;

locDatacopy=obj.locData.copy;
locDatacopy.regroup(150,minframes);

beadind=find(locDatacopy.grouploc.numberInGroup>minframes);
xb=locDatacopy.grouploc.xnm(beadind);
yb=locDatacopy.grouploc.ynm(beadind);
filenumber=locDatacopy.grouploc.filenumber(beadind);
winsize=250;
for k=length(beadind):-1:1
%     beadhere=(mywithin(locDatacopy.loc.xnm,[xb(k)-winsize/2 winsize],locDatacopy.loc.ynm,[yb(k)-winsize/2 winsize]));
%     beadfile=locDatacopy.loc.filenumber==filenumber(k);
%     beadhere=beadhere&beadfile;
    bead(k).loc=locDatacopy.getloc({'xnm','ynm','PSFxnm','PSFynm','frame','phot','znm'},'filenumber',filenumber(k),'Position',[xb(k) yb(k) winsize]);
%     bead(k).loc=copystructReduce(locDatacopy.loc,beadhere,{'xnm','ynm','PSFxnm,','PSFynm','frame'});
    bead(k).filenumber=filenumber(k);
    bead(k).pos=[xb(k) yb(k)];
    fn=fieldnames(bead(k).loc);
    for f=1:length(fn)
        bead(k).loc.(fn{f})=double(bead(k).loc.(fn{f}));
    end    
end
end
