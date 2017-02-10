function imageslicer(varargin)
%optional: x,y,z
%V
%optional handle

%input parser: x,y,z (optional),V,'name',value
%Parent,all gui parameters: scale,  contrastmode, contrast,

%todo: rgb, more dimension, choose dim1,dim2
p=parseinput(varargin);
% extend to rgb
if ~isempty(p.x)
    V=p.z;
    axl{1}=p.V;
    axl{2}=p.x;
    axl{3}=p.y;
else   
    V=p.V;
    axl{1}=1:size(V,1);
    axl{2}=1:size(V,2);
    axl{3}=1:size(V,3);
end

if isempty(p.Parent)
    phandle=figure;
else
    phandle=p.Parent;
end

maxV=max(V(:));

 dim=1:2;
 dimmenu=3;
dimrgb=[];
if isprop(phandle,'KeyPressFcn')
    phandle.KeyPressFcn=@keypress;
end
ax=axes('Parent',phandle,'Position',[0.05,0.18,.95,.82]);
ax.XLim=[0 Inf];

vp1=0.07;
vp2=0.02;
hs=uicontrol('Parent',phandle,'Style','slider','Units','normalized','Position',[0.05 vp1 0.35 0.05],...
    'Min',1,'Max',size(V,3),'Value',1,'SliderStep',[1/(size(V,3)-1) 5/(size(V,3)-1)],'Callback',{@slidercallback,1});
hs2=uicontrol('Parent',phandle,'Style','slider','Units','normalized','Position',[0.05 vp2 0.35 0.05],...
    'Min',1,'Max',size(V,3),'Value',1,'SliderStep',[1/(size(V,3)-1) 5/(size(V,3)-1)],'Callback',{@slidercallback,1});
hn=uicontrol('Parent',phandle,'Style','edit','Units','normalized','String','1','Position',[0.4 vp1 0.075 0.05],'Callback',{@framecallback,1});
hn2=uicontrol('Parent',phandle,'Style','edit','Units','normalized','String','1','Position',[0.4 vp2 0.075 0.05],'Callback',{@framecallback,2});

hmenux=uicontrol('Parent',phandle,'Style','popupmenu','Units','normalized','String',{'x','y','z'},'Position',[0.475 vp1 0.125 0.05],...
    'Callback',@changeaxis);
hmenuy=uicontrol('Parent',phandle,'Style','popupmenu','Units','normalized','String',{'x','y','z'},'Position',[0.6 vp1 0.125 0.05],...
    'Callback',@changeaxis,'Value',2);
haxscale=uicontrol('Parent',phandle,'Style','checkbox','Units','normalized','String','fill','Position',[0.475 vp2 0.1 0.05],...
    'Callback',@plotimage,'Value',1);
hlut=uicontrol('Parent',phandle,'Style','popupmenu','Units','normalized','String',{'parula','gray','hot','jet'},'Position',[0.725 vp1 0.175 0.05],...
    'Callback',@plotimage);
hcontrastcheck=uicontrol('Parent',phandle,'Style','checkbox','Units','normalized','String','global contrast','Position',[0.6 vp2 0.2 0.05],...
    'Callback',@plotimage,'Value',0);
hcontrast=uicontrol('Parent',phandle,'Style','edit','Units','normalized','String','1','Position',[0.8 vp2 0.1 0.05],'Callback',@plotimage);

hresetax=uicontrol('Parent',phandle,'Style','pushbutton','Units','normalized','String','reset','Position',[0.9 vp2 0.1 0.05],'Callback',@resetax);

hrgb=uicontrol('Parent',phandle,'Style','checkbox','Units','normalized','String','RGB','Position',[0.9 vp1 0.1 0.05],'Callback',@updateall);
updateall

plotimage


    function resetax(a,b)
        ax.XLim=[0 Inf];
        plotimage;
    end
    function slidercallback(a,b,slider)
        setslice(a.Value,slider);
    end
    function framecallback(a,b,slider)
       setslice(str2double( a.String),slider);
    end

    function changeaxis(a,b,ax)
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
    
    function setslice(frame,slider)
        frame=min(frame,hs.Max);
        frame=max(1,frame);
        hn.String=num2str(round(frame));
        hs.Value=round(frame);
        plotimage
    end
    function updateall(a,b)
        
        s=size(V);
            
%         for k=1:length(s)
%             dimall{k}=1:s(k);
%         end
        
        if hrgb.Value&&any(s==3)
            dimrgb=find(s==3,1,'last');
            dims=setdiff(3:length(s),dimrgb); 
        else
            dimrgb=[];
            dims=3:length(s);
        end
        str={'x','y'};
        for k=1:length(dims)
            str{end+1}=num2str(dims(k));
        end
        hmenux.String=str;
        hmenuy.String=str;
        
        dim(1)=hmenux.Value;
            dim(2)=hmenuy.Value; 
        dimmenu=setdiff(1:length(s),[dim dimrgb]);
            
        hs.SliderStep=[1/(size(V,dimmenu(1))-1) 5/(size(V,dimmenu(1))-1)];
        hs.Max=size(V,dimmenu(1));
        if length(dimmenu)>1
            hs2.SliderStep=[1/(size(V,dimmenu(2))-1) 5/(size(V,dimmenu(2))-1)];
            hs2.Max=size(V,dimmenu(2));
            setslice(min(hs.Max,hs.Value));
            hs2.Visible='on';
        else
            hs2.Visible='off';
        end
        
    end
    function plotimage(a,b,c)
         s=size(V);
            
        for k=1:length(s)
            dimall{k}=1:s(k);
        end
        xlimold=ax.XLim;
        ylimold=ax.YLim;
        
%         slicez=round(str2double(hn.String));
%         slicet=round(str2double(hn2.String));
        
        if dim(1)<=3
            a1=axl{dim(1)};
        else
            a1=1:size(V,dim(1));
        end
        if dim(2)<=3
            a2=axl{dim(2)};
        else
            a1=1:size(V,dim(2));
        end
        dimall{dimmenu(1)}=round(str2double(hn.String));
        if length(dimmenu)>1
            dimall{dimmenu(2)}=round(str2double(hn2.String));
        end
        img=permute(squeeze(V(dimall{:})),[dim(1) dim(2) dimrgb]);
        
%         switch hmenux.String{hmenux.Value}
%             case 'xy'
%                img=V(:,:,slice)';
%                a1=x;a2=y;
%             case 'yz'
%                img=squeeze(V(slice,:,:))';
%                a1=y;a2=z;
%             case 'xz'
%                img=squeeze(V(:,slice,:))';
%                a1=x;a2=z;
%         end
        
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


function pv=parseinput(in)
p=inputParser;
p.addOptional('x',[],@isnumeric);
p.addOptional('y',[],@isnumeric);
p.addOptional('z',[],@isnumeric);
p.addRequired('V',@isnumeric);
p.addParameter('Parent',[]);
p.addParameter('fill',true,@islogical);
p.addParameter('view',[],@(x) any(validatestring(x,{'xy','xz','yz'})));
parse(p,in{:});
pv=p.Results;

end