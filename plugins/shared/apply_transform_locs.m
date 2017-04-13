function loco=apply_transform_locs(loc,transform,file,p)
if nargin<4||~isfield(p,'datapart') || isempty(p)
    p.datapart.selection='all (T->R)';
end
if p.transformz&&isfield(loc,'znm')
    transformz=true;
else
    transformz=false;
end

if isfield(loc,'xnm_preT')
    loco.xnm_preT=loc.xnm_preT;
    loco.ynm_preT=loc.ynm_preT;
    if transformz
        loco.znm_preT=loc.znm_preT;
    end
else
    loco.xnm_preT=loc.xnm;
    loco.ynm_preT=loc.ynm;
    if transformz
        loco.znm_preT=loc.znm;
    end
end
loco.xnm=loc.xnm;
loco.ynm=loc.ynm;


        
    
loco.frame=loc.frame;
loco.channel=loc.channel;
%only correct file
indf=loc.filenumber==file.number;

if transformz
    loco.znm=loc.znm;
    zf=double(loc.znm);
    indfz=indf;
else
    zf=[];
    indfz=[];
end


xf=double(loc.xnm);
yf=double(loc.ynm);
switch p.datapart.selection
    case {'all (T->R)','all'}
        [x,y,z]=transform.transformCoordinatesInv(xf(indf),yf(indf),zf(indfz));
        indt=true(size(xf));
    case 'all (R->T)'
        [x,y,z]=transform.transformCoordinatesFwd(xf(indf),yf(indf),zf(indfz));
        indt=true(size(xf));    
    case 'target'
        indt=~transform.getRef(xf,yf);
        if ~isempty(indfz)
            indfz=indfz&indt;
        end
        [x,y,z]=transform.transformCoordinatesInv(xf(indf&indt),yf(indf&indt),zf(indfz));
    case 'reference'
        indt=transform.getRef(xf,yf);
        if ~isempty(indfz)
            indfz=indfz&indt;
        end
        [x,y,z]=transform.transformCoordinatesFwd(xf(indf&indt),yf(indf&indt),zf(indfz));
end

indff=find(indf&indt);
loco.xnm(indff)=single(x);
loco.ynm(indff)=single(y);

if transformz
    loco.znm(indff)=single(z);
end

if isfield(p,'setchannel')&&p.setchannel
    loco.ch_preT=loc.channel;
    loco.channel(indf&~indt)=1;
    loco.channel(indff)=2;
end

if isfield(p,'addchannelc')&&p.addchannelc
    loco.ch_preT=loc.channel;
%     loco.channel(indf&~indt)=1;
    loco.channel(indff)=loco.channel(indff)+p.addchannel;
end

