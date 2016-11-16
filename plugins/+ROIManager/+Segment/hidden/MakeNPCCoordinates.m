function l=MakeNPCCoordinates

lpc=2;
dth=pi/48;
theta=0:pi/4:2*pi-eps;
for k=1:lpc-1
    theta=[theta theta+k*dth];
end
[l.x,l.y]=pol2cart(theta,50);
