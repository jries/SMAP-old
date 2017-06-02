function [thetar,phir]=rotateSphericalCoordinates(thetai,phi,axis,alpha)
switch axis
    case 1 %x
        R=[cos(alpha/2) -1i*sin(alpha/2)
            -1i*sin(alpha/2) cos(alpha/2)];
    case 2 %y
         R=[cos(alpha/2) -sin(alpha/2)
            sin(alpha/2) cos(alpha/2)];
        
    case 3 %z
         R=[exp(-1i*alpha/2) 0
            0 exp(1i*alpha/2)];
end


[thetar,phir]=arrayfun(@rotatei,thetai,phi);



function [thetar,phir]=rotatei(thetai,phi)
theta2=pi/4-thetai/2;
z1=[cos(theta2); exp(1i*phi).*sin(theta2)];




z2=R*z1;
thetar=2*atan2(abs(z2(2)),abs(z2(1)));
phir=angle(z2(2))-angle(z2(1));

thetar=pi/2-thetar;
end

end