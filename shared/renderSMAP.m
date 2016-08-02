function imageo=renderSMAP(locs,p,layer,indin,transparency)
%locs: interface.LocalizationData OR locs.xnm/ynm/locprecnm OR locs.x,y,s
%(then take as is)
if nargin==0
    %input parameters
    imageo={'ch_filelist','sr_pixrec','sr_axes','sr_pos','sr_size','rendermode','render_colormode',...
                'renderfield','colorfield_min','colorfield_max','groupcheck','lut','shiftxy_min','shiftxy_max'...
                'mingaussnm','mingausspix','gaussfac','sr_sizeRecPix','shift','displayLocsInd','cam_pixelsize_nm','remout',...
                'rangex','rangey','intensitycoding','sr_layersseparate','sr_layerson'};
    return          
end

if nargin<5
    transparency.mode=1;
end
if strcmpi('tiff', p.rendermode.selection)%obj.locData.files.file(p.ch_filelist.value).istiff
    file=locs.files.file;
    imageo=tif2srimage(file,p);
    imageo.istiff=1;
    return
end

if strcmpi('raw', p.rendermode.selection) %obj.locData.files.file(p.ch_filelist.value).istiff
    file=locs.files.file;
    imageo=tif2srimage(file,p,'raw');
    imageo.istiff=1;
    return
end

if nargin>3&&isempty(indin)
    imageo.image=[];
    imageo.numberOfLocs=0;
    imageo.istiff=0;
    imageo.rangex=[];
    imageo.rangey=[];
    return
end

if nargin<3||isempty(layer)
    layer=p.layer;
end

if length(p.sr_pixrec)==1;
    p.sr_pixrec(2)=p.sr_pixrec(1);
end

imageo.istiff=0;
if isa(locs,'interfaces.LocalizationData')
%        locsh=locs.getloc({'x','y','sx','sy','c','s',p.renderfield.selection},'layer',layer,'position','default');
    locsh=locs.getloc({'xnm','ynm','znm','x','y','locprecnm','sx','sy','PSFxnm','c','s','intensity_render','phot','numberInGroup',p.renderfield.selection},'layer',layer,'position','default');
else
    locsh=locs;
end


if nargin<4||isempty(indin)
    fn=fieldnames(locsh);
    indin=true(length(locsh.(fn{1})),1);
end
    
if isfield(p,'sr_layersseparate')&&~isempty(p.sr_layersseparate)&&p.sr_layersseparate
    if isfield(p,'sr_size')&&~isempty(p.sr_size)&&isfield(p,'sr_layerson')
        p.sr_size(1)=p.sr_size(1)/sum(p.sr_layerson);
    end
end
    
if ~isfield(locsh,'x')||isempty(locsh.x)
%     length(locsh.xnm)
%     length(indin)
    pos.x=locsh.xnm(indin);
else
    pos.x=locsh.x(indin);
end
if ~isfield(locsh,'y')||isempty(locsh.y)
    pos.y=locsh.ynm(indin);
else
    pos.y=locsh.y(indin);
end

if ~isfield(p,'rangex')||isempty(p.rangex)
    rangex=[p.sr_pos(1)-p.sr_size(1) p.sr_pos(1)+p.sr_size(1)];
else
    rangex=p.rangex;
end
if ~isfield(p,'rangey')||isempty(p.rangey)
    rangey=[p.sr_pos(2)-p.sr_size(2) p.sr_pos(2)+p.sr_size(2)];
else   
    rangey=p.rangey;
end
dl=0;

rangexrec=rangex-p.shiftxy_min;
rangeyrec=rangey-p.shiftxy_max; 

if isfield(locsh,'intensity_render')&&~isempty(locsh.intensity_render)
    pos.N=locsh.intensity_render(indin);
elseif isfield(p,'intensitycoding')
    switch p.intensitycoding.selection
        case 'photons'
            pos.N=locsh.phot;
        case 'blinks'
            pos.N=locsh.numberInGroup;
    end
end

lutall=mymakelut(p.lut.selection);
switch p.render_colormode.selection
    case 'normal'
        lut=0;
    case 'z'
        lut=lutall;
        pos.c=locsh.znm(indin);
    case 'field'
        lut=lutall;
        pos.c=locsh.(p.renderfield.selection)(indin);
end

if p.remout==0&&isfield(pos,'c')&&isfield(p,'colorfield_max')&&isfield(p,'colorfield_min')
    pos.c(pos.c>p.colorfield_max)=p.colorfield_max;
    pos.c(pos.c<p.colorfield_min)=p.colorfield_min;
end     
tpar=0;

%build pos
switch lower(p.rendermode.selection)
    case {'gauss','other','constgauss'}
        if ~isfield(locsh,'sx') || ~isfield(locsh,'sy') || isempty(locsh.sx)|| isempty(locsh.sy) 
            if isfield(locsh,'s')&&~isempty(locsh.s)
                sd=locsh.s(indin);
            else
                sd=locsh.locprecnm(indin);
            end
            pos.s=min(max(sd*p.gaussfac,max(p.mingaussnm,p.mingausspix*p.sr_pixrec(1))),400);
        else
            pos.sx=max(locsh.sx(indin)*p.gaussfac,max(p.mingaussnm,p.mingausspix*p.sr_pixrec(1)));
            pos.sy=max(locsh.sy(indin)*p.gaussfac,max(p.mingaussnm,p.mingausspix*p.sr_pixrec(2)));
            switch transparency.mode
                case 1
%             if isempty(transparency)
           
                case 2
     
                    tpar=transparency.parameter;
                case 3
                    if isfield(locsh,'ballradius')
                        pos.sx=locsh.ballradius(indin);
                    else
                      pos.sx(:)=transparency.parameter(2);
                    end
                     pos.sy=pos.sx;
      
                    tpar=transparency.parameter;
            end
           
        end
        if isfield(locsh,'perspective')
            if isfield(pos,'sx')
                pos.sx=pos.sx.*locsh.perspective(indin);
            end
            if isfield(pos,'sy')
                pos.sx=pos.sy.*locsh.perspective(indin);
            end
            if isfield(pos,'s')
                pos.sx=pos.s.*locsh.perspective(indin);
            end
        end
    case 'dl'
        if ~isfield(p,'cam_pixelsize_nm')
            p.cam_pixelsize_nm=140;
        end
        p.sr_pixrec=p.cam_pixelsize_nm;
        if isfield(locsh,'PSFxnm')&&~isempty(locsh.PSFxnm)
            pos.s=locsh.PSFxnm(indin);
        else
            pos.s=pos.x*0+p.cam_pixelsize_nm;
        end
 
        dl=1;
        imageo.istiff=1;
%         rangex(1:2)=rangex(1:2)-p.sr_pixrec/2;
%         rangey(1:2)=rangey(1:2)-p.sr_pixrec/2;
    case {'hist'}
        frender=@histrender;
    case 'tiff'
    

end

% set frender
switch lower(p.rendermode.selection)
    case {'gauss','other'}
        if ~isfield(pos,'sx') 
            frender=@gaussrender;
        else    
            switch transparency.mode
                case 1
                    frender=@gaussrender_ellipt;
                case 2
                    frender=@gaussrenderT_ellipt;
                case 3
                    frender=@gaussrenderTx_ellipt;
            end
           
        end
    case 'dl'
        frender=@gaussrender;
    case {'hist'}
        frender=@histrender;
    case {'constgauss'}
        frender=@constgaussrender;
    case 'tiff'
    

end


rangexrec=rangexrec-p.sr_pixrec(1)/2;
rangeyrec=rangeyrec-p.sr_pixrec(end)/2;

% tic
[srimage,nlocs]=frender(pos,rangexrec, rangeyrec, p.sr_pixrec(1), p.sr_pixrec(end),lut,[p.colorfield_min p.colorfield_max],tpar);
% toc
if dl
    if isfield(p,'sr_sizeRecPix')
        newsize=round(p.sr_sizeRecPix([2 1]));
        srimage=imresize(srimage,[newsize(1) newsize(2)],'nearest');

    end

end

imageo.image=srimage;
imageo.lut=lutall;
imageo.rangex=rangex;
imageo.rangey=rangey;
imageo.numberOfLocs=nlocs;

end
