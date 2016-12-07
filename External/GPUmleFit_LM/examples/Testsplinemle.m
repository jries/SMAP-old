clear all

load('matlab161027.mat')
coeff=single(coeff);

Nfits = 100;
Nphotons = 500;
Npixels = 13;
bg = 10;

coordsxy = Npixels/2 -1 +rand([Nfits 2]);
coordsz = 25*rand(Nfits,1);

output = single(zeros(Npixels,Npixels,Nfits));
spline_xsize = size(coeff,1);
spline_ysize = size(coeff,2);
spline_zsize = size(coeff,3);
off = ((spline_xsize+1)-2*Npixels)/2;


for kk = 1:Nfits
    xcenter = coordsxy(kk,1);
    ycenter = coordsxy(kk,2);
    zcenter = coordsz(kk);
    xc = -2*(xcenter - 6.5+0.5);
    yc = -2*(ycenter - 6.5+0.5);
    zc = zcenter -floor(zcenter);
    xstart = 0;
    while xc>1
        xstart = xstart+1;
        xc = xc-1;
    end
    
    while xc<0
        xstart = xstart-1;
        xc = xc+1;
    end
    
    ystart = 0;
    while yc>1
        ystart = ystart+1;
        yc = yc-1;
    end
    
    while yc<0
        ystart = ystart-1;
        yc = yc+1;
    end
    
    zstart = floor(zcenter);
    [delta_f,delta_dxf,delta_ddxf,delta_dyf,delta_ddyf,delta_dzf,delta_ddzf]=computeDelta3D(single(xc),single(yc),single(zc));
    for ii = 0:Npixels-1
        for jj = 0:Npixels-1
             temp = fAt3D(2*ii+xstart+off,2*jj+ystart+off,zstart,spline_xsize,spline_ysize,spline_zsize,delta_f,coeff);
             model = temp*Nphotons+bg;
             output(ii+1,jj+1,kk)=model;
            
        end
    end

    kk
    
end
outputtemp = repmat(output,1,1,1000);
% outputtemp=single(noise(outputtemp,'poisson',1));
outputtemp = single(poissrnd(outputtemp));
[P_Spline_sim CRLB LL]=GPUmleFit_LM(single(outputtemp(:,:,1:10000)),single(coeff),50,5,0);

resultfit = [P_Spline_sim(:,1:2) P_Spline_sim(:,5)];
resultsim = [coordsxy coordsz];
resultfit(1:100,:)-resultsim