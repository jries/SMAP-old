function imageo=renderSMAP(locs,p,layer,indin,transparency)
%locs: interface.LocalizationData OR locs.xnm/ynm/locprecnm OR locs.x,y,s
%(then take as is)
if nargin==0
    %input parameters
    imageo={'ch_filelist','sr_pixrec','sr_axes','sr_pos','sr_size','rendermode','renderp1',...
                'renderfield','colorfield_min','colorfield_max','groupcheck','lut','shiftxy_min','shiftxy_max'...
                'mingaussnm','mingausspix','gaussfac','sr_sizeRecPix','shift','displayLocsInd','cam_pixelsize_nm','remout',...
                'rangex','rangey'};
    return          
end
if nargin<5
    transparency=[];
end
if strcmpi('tiff', p.rendermode.selection) %obj.locData.files.file(p.ch_filelist.value).istiff
    file=locs.files.file;
    imageo=tif2srimage(file,p);
    imageo.istiff=1;
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
    locsh=locs.getloc({'xnm','ynm','znm','x','y','locprecnm','sx','sy','PSFxnm','c','s','intensity_render',p.renderfield.selection},'layer',layer,'position','default');
else
    locsh=locs;
end


if nargin<4
    fn=fieldnames(locsh);
    indin=true(length(locsh.(fn{1})),1);
end
    
    
if isempty(locsh.x)
    pos.x=locsh.xnm(indin);
else
    pos.x=locsh.x(indin);
end
if isempty(locsh.y)
    pos.y=locsh.ynm(indin);
else
    pos.y=locsh.y(indin);
end

if ~isfield(p,'rangex')||isempty(p.rangex)
    rangex=[p.sr_pos(1)-p.sr_size(1) p.sr_pos(1)+p.sr_size(1)]-p.shiftxy_min;
else
    rangex=p.rangex;
end
if ~isfield(p,'rangey')||isempty(p.rangey)
    rangey=[p.sr_pos(2)-p.sr_size(2) p.sr_pos(2)+p.sr_size(2)]-p.shiftxy_max; 
else   
    rangey=p.rangey;
end
dl=0;

if isfield(locsh,'intensity_render')&&~isempty(locsh.intensity_render)
    pos.N=locsh.intensity_render(indin);
end

lutall=mymakelut(p.lut.selection);
switch p.renderp1.selection
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

switch lower(p.rendermode.selection)
    case {'gauss','other'}
        if isempty(locsh.sx)|| isempty(locsh.sy) 
            if ~isempty(locsh.s)
                sd=locsh.s(indin);
            else
                sd=locsh.locprecnm(indin);
            end
            pos.s=min(max(sd*p.gaussfac,max(p.mingaussnm,p.mingausspix*p.sr_pixrec(1))),400);
            frender=@gaussrender;
        else
            pos.sx=max(locsh.sx(indin)*p.gaussfac,max(p.mingaussnm,p.mingausspix*p.sr_pixrec(1)));
            pos.sy=max(locsh.sy(indin)*p.gaussfac,max(p.mingaussnm,p.mingausspix*p.sr_pixrec(2)));
            if isempty(transparency)
                frender=@gaussrender_ellipt;
             
            else
                frender=@gaussrenderT_ellipt;
            end
           
        end
    case 'dl'
        if ~isfield(p,'cam_pixelsize_nm')
            p.cam_pixelsize_nm=140;
        end
        p.sr_pixrec=p.cam_pixelsize_nm;
        pos.s=loch.PSFxnm(indin);
        if isempty(pos.s)
            pos.s=pos.x*0+p.cam_pixelsize_nm;
        end
        frender=@gaussrender;
        dl=1;
    case 'hist'
        frender=@histrender;
    case 'tiff'
    

end
            
          


[srimage,nlocs]=frender(pos,rangex, rangey, p.sr_pixrec(1), p.sr_pixrec(end),lut,[p.colorfield_min p.colorfield_max],transparency);
if dl
    if isfield(p,'sr_sizeRecPix')
        newsize=round(p.sr_sizeRecPix([2 1]));
        srimage=imresize(srimage,[newsize(1) newsize(2)],'nearest');

    end

end

imageo.image=srimage;
imageo.lut=lutall;
imageo.rangex=rangex+p.shiftxy_min;
imageo.rangey=rangey+p.shiftxy_max;
imageo.numberOfLocs=nlocs;

end
