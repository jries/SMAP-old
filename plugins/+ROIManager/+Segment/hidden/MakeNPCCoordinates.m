function l=MakeNPCCoordinates
radius=50;
lpc=4; %number of corners
dth=pi/128; %shift between corners
theta=0:pi/4:2*pi-eps; %8 corners
thetaa=[];
for k=1:lpc-1
    thetaa=[thetaa theta+k*dth];
end
[l.x,l.y]=pol2cart(thetaa,radius);
