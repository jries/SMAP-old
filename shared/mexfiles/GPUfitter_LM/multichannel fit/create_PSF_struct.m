PSF = struct('ipalm_im', nan, 'I', nan, 'A', nan, 'B', nan, 'Ispline', nan, ...
    'Aspline', nan, 'Bspline', nan, 'k', nan, 'zstep', nan, 'zcand', nan);


%load('ipalm_im_no_ab_x'); %load zstack
load('ipalm_im_no_ab_x');

ipalm_im = ipalm_im * 4 / sum(reshape(ipalm_im(:, :, round(length(zcand) / 2), :), 1, [])); %normalize PSF
imageslicer(ipalm_im);

%phaseshift, lambdanm and zcand (z range in nm) are 
%saved in the file with the zstack
PSF.phaseshift = phaseshift;
PSF.k = 2 * pi / lambdanm;
PSF.zcand = zcand;
PSF.zstep = zcand(2) - zcand(1);

[x_size, y_size, z_size, phi_size] = size(ipalm_im);
PSF.ipalm_im = [];
sum_num = 4; % number of neighbouring pixels to sum up
if sum_num == 1
    PSF.ipalm_im = ipalm_im;
else
    for ii = 1 : 1 : floor(x_size / sum_num)
        for jj = 1 : 1 : floor(y_size / sum_num)
            tmp = sum(ipalm_im(ii * sum_num - sum_num + 1 : ii * sum_num, jj * sum_num - sum_num + 1 : jj * sum_num, :, :));
            PSF.ipalm_im(ii, jj, :, :) = sum(tmp);
        end
    end
end
imageslicer(squeeze(PSF.ipalm_im(:, :, :, 1)));
imageslicer(squeeze(PSF.ipalm_im(:, :, :, 2)));
imageslicer(squeeze(PSF.ipalm_im(:, :, :, 3)));
imageslicer(squeeze(PSF.ipalm_im(:, :, :, 4)));
subsz = floor(x_size / sum_num);

%make I, A, B and their splines
index_0 = find(PSF.phaseshift == 0);
index_pi = find(PSF.phaseshift == pi);
index_pi2 = find(PSF.phaseshift == pi/2);

PSF.I = squeeze((PSF.ipalm_im(:, :, :, index_0) + PSF.ipalm_im(:, :, :, index_pi)) / 2);

kz2 = 2 * PSF.k * PSF.zcand';
kz2 = permute(repmat(kz2, 1, subsz, subsz), [2, 3, 1]);
cos2kz = cos(kz2);
sin2kz = sin(kz2);

F_0 = squeeze(PSF.ipalm_im(:, :, :, index_0)) - PSF.I;
F_pi2 = squeeze(PSF.ipalm_im(:, :, :, index_pi2)) - PSF.I;

PSF.A = F_0 .* cos2kz - F_pi2 .* sin2kz;
PSF.B = F_0 .* sin2kz + F_pi2 .* cos2kz;
%imageslicer(PSF.A);
%imageslicer(PSF.B);

PSF.Ispline = Spline3D_interp(PSF.I);
PSF.Aspline = Spline3D_interp(PSF.A);
PSF.Bspline = Spline3D_interp(PSF.B);

save('PSF');


