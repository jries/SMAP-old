function lut=mymakelut(lutname)
if nargin==0
    lut={'red hot','green cold','cyan cold', 'red','green','blue','cyan','jet','bgy','bry','gray','gray inv','hsv','parula'};
else
switch lutname
    case 'red hot'
        lut=hot(256);
    case 'red'
        lut=zeros(256,3);
        lut(:,1)=(0:255)/255;
    case 'green'
        lut=zeros(256,3);
        lut(:,2)=(0:255)/255;
    case 'blue'
        lut=zeros(256,3);
        lut(:,3)=(0:255)/255;    
    case 'cyan'
        lut=zeros(256,3);
        lut(:,3)=(0:255)/255; 
        lut(:,2)=(0:255)/255; 
    case 'jet'
        lut=jet(256);
    case 'green cold'
        lut=usercolormap([0 0 0],[0 1 0],[0 1 1],[1 1 1]);
    case 'cyan cold'
        lut=usercolormap([0 0 0],[0 .2 1],[0 0.5 1],[0 1 1],[.5 1 1],[1 1 1]);
    case 'bgy'
        lut=usercolormap([0 0 1],[0 1 0],[1 1 0]);
    case 'bry'
        lut=usercolormap([0 0 1],[1 0 0],[1 1 0]);
    case 'gray'
        lut=gray(256);
    case 'gray inv'
        lut=gray(256);
        lut=lut(end:-1:1,:);
    case 'hsv'
        lut=hsv(256);
    case 'parula'
        lut=parula(256);

end
end