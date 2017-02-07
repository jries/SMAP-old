function imageslicer(varargin)
%optional: x,y,z
%V
%optional handle
if nargin>2
    V=varargin{4};
    x=varargin{1};
    y=varargin{2};
    z=varargin{3};
    if nargin>4
        handle=varargin{5};
    else
        handle=figure;
    end
else   
    V=varargin{1};
    x=1:size(V,1);
    y=1:size(V,2);
    z=1:size(V,3);
    if nargin>1
        handle=varargin{2};
    else
        handle=figure;
    end  
end

maxV=max(V(:));

dim=1:3;

if isprop(handle,'KeyPressFcn')
handle.KeyPressFcn=@keypress;
end
ax=axes('Parent',handle,'Position',[0.05,0.18,.95,.82]);
ax.XLim=[0 Inf];

vp1=0.07;
vp2=0.02;
hs=uicontrol('Parent',handle,'Style','slider','Units','normalized','Position',[0.05 vp1 0.35 0.05],...
    'Min',1,'Max',size(V,3),'Value',1,'SliderStep',[1/(size(V,3)-1) 5/(size(V,3)-1)],'Callback',@slidercallback);
hn=uicontrol('Parent',handle,'Style','edit','Units','normalized','String','1','Position',[0.4 vp1 0.1 0.05],'Callback',@framecallback);
hmenu=uicontrol('Parent',handle,'Style','popupmenu','Units','normalized','String',{'xy','yz','xz'},'Position',[0.05 vp2 0.15 0.05],...
    'Callback',@changeaxis);
haxscale=uicontrol('Parent',handle,'Style','checkbox','Units','normalized','String','fill','Position',[0.2 vp2 0.1 0.05],...
    'Callback',@plotimage,'Value',1);
hlut=uicontrol('Parent',handle,'Style','popupmenu','Units','normalized','String',{'parula','gray','hot','jet'},'Position',[0.725 vp1 0.175 0.05],...
    'Callback',@plotimage);
hcontrastcheck=uicontrol('Parent',handle,'Style','checkbox','Units','normalized','String','global contrast','Position',[0.6 vp2 0.2 0.05],...
    'Callback',@plotimage,'Value',0);
hcontrast=uicontrol('Parent',handle,'Style','edit','Units','normalized','String','1','Position',[0.8 vp2 0.1 0.05],'Callback',@plotimage);

hresetax=uicontrol('Parent',handle,'Style','pushbutton','Units','normalized','String','reset','Position',[0.3 vp2 0.1 0.05],'Callback',@resetax);


plotimage
    function resetax(a,b)
        ax.XLim=[0 Inf];
        plotimage;
    end
    function slidercallback(a,b)
        setslice(a.Value);
    end
    function framecallback(a,b)
       setslice(str2double( a.String));
    end

    function changeaxis(a,b)
        switch a.String{a.Value}
            case 'xy'
                dim = [1 2 3];
            case 'yz'
                dim = [2 3 1];
            case 'xz'
                dim = [1 3 2];
        end
        ax.XLim=[0 Inf];
        updateall

    end
    
    function setslice(frame)
        frame=min(frame,hs.Max);
        frame=max(1,frame);
        hn.String=num2str(round(frame));
        hs.Value=round(frame);
        plotimage
    end
    function updateall
        hs.SliderStep=[1/(size(V,dim(3))-1) 5/(size(V,dim(3))-1)];
        hs.Max=size(V,dim(3));
        setslice(min(hs.Max,hs.Value));
        
    end
    function plotimage(a,b,c)
        
        xlimold=ax.XLim;
        ylimold=ax.YLim;
        slice=round(str2double(hn.String));
        switch hmenu.String{hmenu.Value}
            case 'xy'
               img=V(:,:,slice)';
               a1=x;a2=y;
            case 'yz'
               img=squeeze(V(slice,:,:))';
               a1=y;a2=z;
            case 'xz'
               img=squeeze(V(:,slice,:))';
               a1=x;a2=z;
        end
        
        %contrast
        if hcontrastcheck.Value
            imax=str2double(hcontrast.String)*maxV;
        else
            imax=str2double(hcontrast.String)*max(img(:));
        end
        
        img(img>imax)=imax;
        img(1,1)=imax;
        imagesc(ax,a1,a2,img);
        if haxscale.Value
            axis(ax,'fill')
        else
            axis(ax,'equal')
        end
        if ~isinf(xlimold(2))
        ax.XLim=xlimold;ax.YLim=ylimold;
        end
        colormap(ax,hlut.String{hlut.Value})
        
    end
    function keypress(a,b)
        if strcmp(b.Character,'+')||strcmp(b.Key,'rightarrow')
            frame=hs.Value;
            setslice(frame+1);
        elseif strcmp(b.Character,'-')||strcmp(b.Key,'leftarrow')
            frame=hs.Value;
            setslice(frame-1);
        elseif strcmp(b.Key,'uparrow')
            hcontrast.String=num2str(str2double(hcontrast.String)*1.1,'%1.2f');
            plotimage;
        elseif strcmp(b.Key,'downarrow')
            hcontrast.String=num2str(str2double(hcontrast.String)*.9,'%1.2f');
            plotimage;
        end
            
        
    end
end