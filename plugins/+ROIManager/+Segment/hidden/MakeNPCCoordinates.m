theta=0:pi/4:2*pi-eps;
[x,y]=pol2cart(theta,50);
[f,p]=uiputfile('*.mat');
save([p f],'x','y')