function loco=apply_transform_locs(loc,transform,file,p)
if nargin<4||~isfield(p,'datapart') || isempty(p)
    p.datapart.selection='all';
end
loco.xnm_preT=loc.xnm;
loco.ynm_preT=loc.ynm;
loco.xnm=loc.xnm;
loco.ynm=loc.ynm;
loco.frame=loc.frame;
%only correct file
indf=loc.filenumber==file.number;

xf=double(loc.xnm);
yf=double(loc.ynm);
switch p.datapart.selection
    case 'all'
    [x,y]=transform.transformCoordinatesFwd(xf,yf);
    indt=true(size(x));
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
    loco.channel=ones(size(loc.xnm))*1;
    loco.channel(indff)=2;
end

