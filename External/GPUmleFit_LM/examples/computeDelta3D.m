function [delta_f,delta_dxf,delta_ddxf,delta_dyf,delta_ddyf,delta_dzf,delta_ddzf] = computeDelta3D(z_delta, y_delta, x_delta)

% if x_delta<0
%     x_delta = 0;
% end
% 
% if x_delta>1
%     x_delta = 1;
% end
% 
% if y_delta<0
%     y_delta = 0;
% end
% 
% if y_delta>1
%     y_delta = 1;
% end
% 
% if z_delta<0
%     z_delta = 0;
% end
% 
% if z_delta>1
%     z_delta = 1;
% end

delta_f = single(zeros(64,1));
delta_dxf = single(zeros(64,1));
delta_dyf = single(zeros(64,1));
delta_dzf = single(zeros(64,1));
delta_ddxf = single(zeros(64,1));
delta_ddyf = single(zeros(64,1));
delta_ddzf = single(zeros(64,1));
cx = single(1);
for i = 0:3
    cy = single(1);
    for j = 0:3
        cz = single(1);
        for k = 0:3
            delta_f(i*16+j*4+k+1) = single(cx*cy*cz);
            if k<3
                delta_dxf(i*16+j*4+k+1+1) = single((k+1))*cx*cy*cz;
            end
            if k<2
                delta_ddxf(i*16+j*4+k+2+1)=single((k+1)*(k+2))*cx*cy*cz;
            end
            if j<3
                delta_dyf(i*16+(j+1)*4+k+1) = single((j+1))*cx*cy*cz;
            end
            if j<2
                delta_ddyf(i*16+(j+2)*4+k+1)= single((j+1)*(j+2))*cx*cy*cz;
            end
            if i<3
                delta_dzf((i+1)*16+j*4+k+1) = single((i+1))*cx*cy*cz;
            end
            if i<2
                delta_ddzf((i+2)*16+j*4+k+1)=single((i+1)*(i+2))*cx*cy*cz;
            end
            cz = cz*z_delta;
        end
        cy = cy*y_delta;
    end
    cx = cx*x_delta;
end