function transform=transform_locs(locData,p)

%get fields
if p.uselayers
    [locref,indref]=locData.getloc({'xnm','ynm','frame','filenumber','znm','PSFxnm'},'layer',p.reflayer.Value,'position','roi');
    [loctarget,indtarget]=locData.getloc({'xnm','ynm','frame','filenumber','znm','PSFxnm'},'layer',p.targetlayer.Value,'position','roi');
else
    locref=locData.getloc({'xnm','ynm','frame','filenumber','znm','PSFxnm'});
%     loctarget=locData.getloc('xnm','ynm','frame','filenumber','znm','PSFxnm');
    
    filter=locData.layer(1).filter;
    ind=true(length(locref.xnm),1);
    if isfield(filter,'locprecnm')
        ind=ind&filter.locprecnm;
    end
    if isfield(filter,'PSFxnm')
        ind=ind&filter.locprecnm;
    end
    if isfield(filter,'frame')
        ind=ind&filter.locprecnm;
    end
    ind=ind&locref.filenumber==p.dataselect.Value;
%     indref=ind;
%     indtarget=ind;
    fn=fieldnames(locref);
    for k=1:length(fn)
        if ~isempty(locref.(fn{k}))
        locref.(fn{k})=locref.(fn{k})(ind);
        end
%         loctarget.(fn{k})=loctarget.(fn{k})(ind);
    end
    loctarget=locref;
end

locref.x=locref.xnm;
locref.y=locref.ynm;

loctT=loctarget;
% file=locData.files.file(p.layer1_.ch_filelist.Value);

chipsizenm=p.currentfileinfo.pixsize*512*1000; 
facsize=ones(2,1);
separator=[chipsizenm chipsizenm];
switch p.targetpos.selection
    case 'top'
        dy=-chipsizenm/2;
        dx=0;
        separator(2)=chipsizenm/2;
        indtarget=loctarget.ynm<separator(2);
    case 'bottom'
        dy=chipsizenm/2;
        dx=0;
        separator(2)=chipsizenm/2;
        indtarget=loctarget.ynm>separator(2);
    case 'left'
        dx=-chipsizenm/2;
        dy=0;
        separator(1)=chipsizenm/2;
        indtarget=loctarget.xnm<separator(1);
    case 'right'
        dx=chipsizenm/2;
        dy=0;
        separator(1)=chipsizenm/2;
        indtarget=loctarget.xnm>separator(1);
    otherwise
        dx=0;
        dy=0;
        indtarget=true(size(loctarget.xnm));
end

if p.useT
    Tload=load(p.Tfile);
    Tinitial=Tload.transformation;
    [loctT.x,loctT.y]=Tinitial.transformCoordinatesInv(loctarget.xnm,loctarget.ynm);
    mirrorinfo=Tinitial.tinfo.mirror;
    %     pos=Tinitial.pos;
%     size=Tinitial.size;
else %all initial estimation:
%     approximate shift from size and position
   
    loctT.x=loctT.xnm-dx;
    loctT.y=loctT.ynm-dy;
    %  initial mirror
    midmirror=chipsizenm/4;
    mirrorinfo.midmirror=midmirror;
    mirrorinfo.targetmirror=p.targetmirror.selection;
    
    switch p.targetmirror.selection
        case 'left-right'
            loctT.x=2*midmirror-loctT.x;
            
        case 'up-down' 
            loctT.y=2*midmirror-loctT.y;
            
        case 'both'  
            loctT.x=2*midmirror-loctT.x;
            loctT.y=2*midmirror-loctT.y;
           
    end

end



% indref=loctarget.xnm<chipsizenm*facsize(1)&loctarget.ynm<chipsizenm*facsize(2);
% loctT=copystructReduce(loctT,ind);
loctT.frame(~indtarget)=-1;

pixelsizerec=p.register_parameters.pixelsizenm;
roi=p.currentfileinfo.roi;
roinm=roi*p.currentfileinfo.pixsize*1000;

% pos=[mean(roinm([1,3])) meannm(roi([2,4]))];
%     rsize=[roi(3)-roi(1) roi(4)-roi(2)];

rangex=[roinm(1) roinm(1)+roinm(3)*facsize(1)];
rangey=[roinm(2) roinm(2)+roinm(4)*facsize(2)];


imr=myhist2(locref.x(~indtarget),locref.y(~indtarget),pixelsizerec,pixelsizerec,rangex,rangey);
imt=myhist2(loctT.x(indtarget),loctT.y(indtarget),pixelsizerec,pixelsizerec,rangex,rangey);
% 
% figure(99);
% hold off
% plot(locref.x,locref.y,'bo')
% hold on
% plot(loctT.x,loctT.y,'rx')
% % plot(loctarget.xnm,loctarget.ynm,'g.')
% 
% 


maxshift=round(p.register_parameters.maxshift_corr/pixelsizerec);
ax1=initaxis(p.resultstabgroup,'shiftcorr');

% s=size(imr)
% ima=zeros(s(1),s(2),3);
% ima(:,:,1)=imr;ima(:,:,2)=imt;
% figure(88);imagesc(ima)
% asdf
[dxpt,dypt]=getShiftCorr(imr,imt,1,maxshift);
dxcorr=dxpt*pixelsizerec;
dycorr=dypt*pixelsizerec;

loctT.x=loctT.x+dxcorr;
loctT.y=loctT.y+dycorr;
% [loctarget,indt]=transform.getPartLocs(loc,loc.x(ig),loc.y(ig),'target');
% [locref,indr]=transform.getPartLocs(loc,loc.x(ig),loc.y(ig),'reference');
% figure(99)
% plot(loctT.x(indtarget),loctT.y(indtarget),'gx',locref.x,locref.y,'bo')
% hold off

% locrefred=reduceFieldLength(locref,ig);
% loctargetred=reduceFieldLength(loctarget,ig);

[iAa,iBa,na,nb,nseen]=matchlocsall(locref,loctT,0,0,p.register_parameters.maxshift_match,p.register_parameters.maxlocsused);

% transform.findTransform(locref.x(iAa),locref.y(iAa),loctarget.x(iBa),loctarget.y(iBa))
transform=interfaces.LocTransform;
t.type=p.transform.selection;
t.parameter=p.transformparam;
transform.findTransform(locref.x(iAa),locref.y(iAa),locref.x(iBa),locref.y(iBa),t)
% transform.findTransform(locref.x(iAa),locref.y(iAa),loctT.x(iBa),loctT.y(iBa),t)
% if p.showresults
    initaxis(p.resultstabgroup,'scatter')
    [xa, ya]=transform.transformCoordinatesInv((locref.x(iBa)),(locref.y(iBa)));
%      [xa, ya]=transform.transformCoordinatesInv((loctT.x(iBa)),(loctT.y(iBa)));
%  [xa, ya]=transform.transformCoordinates((loc.x(indr(iBa))),(loc.y(indt(iBa))),'target');
 
   dx=xa-locref.x(iAa);
   dy=ya-locref.y(iAa);
   
%    dx=loctarget.x(iBa)-locref.x(iAa);
%    dx=loctarget.y(iBa)-locref.y(iAa);
%    figure(88)
   dscatter(dx,dy)
   title({['number of anchor points: ' num2str(length(iBa)) ' of ' num2str(nseen)],['dx= ' num2str(std(dx),3) ' nm, dy= ' num2str(std(dy),3) ' nm']});
 ax3=initaxis(p.resultstabgroup,'hist');
 hist(dx,50)


transform.tinfo.targetpos=p.targetpos.selection;
transform.tinfo.separator=separator;
transform.tinfo.mirror=mirrorinfo;

% [f,path]=uiputfile(p.Tfile);
% if f
%output: transform object
%     transform.T=mytform;
%     transform.midpoint=midp;
%     transform.refpart=p.refpart.value;
%     transform.refmirror=0; %%not implemented
%     transform.targetpart=p.targetpart.value;
%     transform.targetmirror=0; %%not implemented
%     save([path filesep f],'transform')
% end
% 
% x=loc.x(ig);
% y=loc.y(ig);
% f=loc.frame(ig);
% fn=loc.filenumber(ig);


% mx=myquantile(x,0.9


% 
% 
% function outim=getsrimage(l,pixrec,m)
% 
% [outim,dxo,dyo]=myhist2(l.x,l.y,pixrec,pixrec,m.x,m.y);
%    
% 
