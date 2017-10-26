%  Copyright (c)2017 Ries Lab, European Molecular Biology Laboratory,
%  Heidelberg.
%  
%  This file is part of GPUmleFit_LM Fitter.
%  
%  GPUmleFit_LM Fitter is free software: you can redistribute it and/or modify
%  it under the terms of the GNU General Public License as published by
%  the Free Software Foundation, either version 3 of the License, or
%  (at your option) any later version.
%  
%  GPUmleFit_LM Fitter is distributed in the hope that it will be useful,
%  but WITHOUT ANY WARRANTY; without even the implied warranty of
%  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%  GNU General Public License for more details.
%  
%  You should have received a copy of the GNU General Public License
%  along with GPUmleFit_LM Fitter.  If not, see <http://www.gnu.org/licenses/>.
%  
%  
%  Additional permission under GNU GPL version 3 section 7
%  
%  If you modify this Program, or any covered work, by
%  linking or combining it with libraries required for interaction
%  with analysis programs such as Igor Pro or Matlab,
%  the licensors of this Program grant you additional permission
%  to convey the resulting work.
%%
function [out] = simSplinePSF(Npixels,PSF,I,bg,cor, phi0)
t=tic;
if (nargin <5)
   error('Minimal usage: simSplinePSF(Npixels,coeff,I,bg,cor)');
end

% if (bg == 0)
%     bg = 10^-10;
% end    

Nfits = size(cor,1);
spline_xsize = size(PSF.Ispline,1);
spline_ysize = size(PSF.Ispline,2);
spline_zsize = size(PSF.Ispline,3);
off = floor(((spline_xsize+1)-Npixels)/2);
data = zeros(Npixels,Npixels,Nfits,'single');


for kk = 1:Nfits
    xcenter = cor(kk,1);
    ycenter = cor(kk,2);
    zcenter = cor(kk,3);
    
    xc = -1*(xcenter - Npixels/2+0.5);
    yc = -1*(ycenter - Npixels/2+0.5);
    zc = zcenter - floor(zcenter);
    
    xstart = floor(xc);
    xc = xc - xstart;
    
    ystart = floor(yc);
    yc = yc - ystart;
    

    zstart = floor(zcenter);
    
   
    [delta_f,delta_dxf,delta_ddxf,delta_dyf,delta_ddyf,delta_dzf,delta_ddzf]=computeDelta3Dj_v2(single(xc),single(yc),single(zc));
    
    for ii = 0:Npixels-1
        for jj = 0:Npixels-1
             phi = cor(kk, 4);
             for p = 1:1:length(phi0)
                 temp = fAt3Dj_4Pi(ii+xstart+off,jj+ystart+off,zstart,spline_xsize,spline_ysize,spline_zsize,delta_f,PSF, phi, phi0(p));
                 model = temp * I + bg(kk, p);
                 data(ii+1, jj+1, kk, p) = model;
                 if temp < 0
                     disp('Neg');
                 end
             end
        end
    end
    if toc(t)>1
        disp(kk/Nfits)
        t=tic;
    end
    
end
for p = 1:1:length(phi0)
    out(:, :, :, p) = (poissrnd(data(:, :, :, p),Npixels,Npixels,Nfits)); 
   %out = data;
end
