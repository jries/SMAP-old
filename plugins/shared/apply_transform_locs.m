function loco=apply_transform_locs(loc,transform,file,p)
if nargin<4||~isfield(p,'datapart') || isempty(p)
    p.datapart.selection='all (T->R)';
end
if isfield(loc,'xnm_preT')
    loco.xnm_preT=loc.xnm_preT;
    loco.ynm_preT=loc.ynm_preT;
else
    loco.xnm_preT=loc.xnm;
    loco.ynm_preT=loc.ynm;
end
loco.xnm=loc.xnm;
loco.ynm=loc.ynm;
loco.frame=loc.frame;
loco.channel=loc.channel;
%only correct file
indf=loc.filenumber==file.number;

xf=double(loc.xnm);
yf=double(loc.ynm);
switch p.datapart.selection
    case {'all (T->R)','all'}
        [x,y]=transform.transformCoordinatesInv(xf(indf),yf(indf));
        indt=true(size(xf));
    case 'all (R->T)'
        [x,y]=transform.transformCoordinatesFwd(xf(indf),yf(indf));
        indt=true(size(xf));    
    case 'target'
        indt=~transform.getRef(xf,yf);
        [x,y]=transform.transformCoordinatesInv(xf(indf&indt),yf(indf&indt));
    case 'reference'
        indt=transform.getRef(xf,yf);
        [x,y]=transform.transformCoordinatesFwd(xf(indf&indt),yf(indf&indt));
end

indff=find(indf&indt);
loco.xnm(indff)=single(x);
loco.ynm(indff)=single(y);


if isfield(p,'setchannel')&&p.setchannel
    loco.ch_preT=loc.channel;
    loco.channel(indf&~indt)=1;
    loco.channel(indff)=2;
end

