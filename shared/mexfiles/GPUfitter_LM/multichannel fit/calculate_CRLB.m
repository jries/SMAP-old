function [CRLB] =  calculate_CRLB(Nfits, PSF, sz, phi0, Theta)
pi = single(3.141592);
spline_xsize = size(PSF.Ispline, 1);
spline_ysize = size(PSF.Ispline, 2);
spline_zsize = size(PSF.Ispline, 3);

PSFSigma = single(1.5);
xc = single(0);
yc = single(0);
zc = single(0);
NV = 6; %number of parameters to optimize

dudt = single(zeros(sz,sz,NV));

M = zeros(NV, NV, Nfits);
Minv = zeros(NV, NV, Nfits);

for tx = 1:Nfits
    xc = single(-1*(Theta(tx, 1) - sz/2+0.5));
    yc = single(-1*(Theta(tx, 2) - sz/2+0.5));
    zc = single(Theta(tx, 5)-floor(Theta(tx, 5)));
    zstart = single(floor(Theta(tx, 5)));
    
    xstart = floor(xc);
    xc = xc-floor(xc);
    
    ystart = floor(yc);
    yc = yc-floor(yc);
    off = floor(((spline_xsize+1)-sz)/2);
    [delta_f,delta_dxf,delta_ddxf,delta_dyf,delta_ddyf,delta_dzf,delta_ddzf]=computeDelta3Dj_v2(single(xc),single(yc),single(zc));
    
    for ii = 1:sz
        for jj = 1:sz
            
            
            [dudt, model] =  kernel_DerivativeSpline_v2(ii+xstart+off, jj+ystart+off, zstart, spline_xsize, spline_ysize, spline_zsize, delta_f, delta_dxf, delta_dyf, delta_dzf, PSF, squeeze(Theta(tx, :)), NV, phi0);
            for kk = 1:1:length(phi0)
                M(:, :, tx) = M(:, :, tx) + dudt(:, kk) * dudt(:, kk)' / model(kk);
            end
        end
    end
    Minv(:, :, tx) = inv(M(:, :, tx));
    for kk = length(squeeze(Theta(tx, :))):-1:1
        CRLB(tx,kk) = Minv(kk,kk,tx);
    end
end

