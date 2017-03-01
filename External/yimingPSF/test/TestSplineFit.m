% clear all
% 
% load('matlab161027.mat')

Nfits = 10000;
% Nphotons = 20297;
% Nphotons =Nphotons*4;
Npixels = 13;
bg = 1;

% coordsxy = Npixels/2 -1 +rand([Nfits 2]);
coordsxy = Npixels/2 -1 +0.5*ones([35 2]);
% coordsz = 25*rand(Nfits,1);
coordsz = linspace(0,25,35);
coordsz = coordsz';

output = single(zeros(Npixels,Npixels,Nfits));
spline_xsize = size(coeff,1);
spline_ysize = size(coeff,2);
spline_zsize = size(coeff,3);
off = ((spline_xsize+1)-2*Npixels)/2;

Nphotons = 1000;
Nphotons =Nphotons*4;

% coordsxy = [6.5,6,5];
% coordsz =5;
for kk = 1:35
    xcenter = coordsxy(kk,1);
    ycenter = coordsxy(kk,2);
    zcenter = coordsz(kk);
    xc = -2*(xcenter - 6.5-0.5);
    yc = -2*(ycenter - 6.5-0.5);
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
%      [delta_f,delta_dxf,delta_ddxf,delta_dyf,delta_ddyf,delta_dzf,delta_ddzf]=computeDelta3D(single(zc),single(yc),single(xc));
%     [delta_f,delta_dxf,delta_ddxf,delta_dyf,delta_ddyf,delta_dzf,delta_ddzf]=computeDelta3D(single(xc),single(yc),single(zc));%good
    [delta_f,delta_dxf,delta_ddxf,delta_dyf,delta_ddyf,delta_dzf,delta_ddzf]=computeDelta3D(single(yc),single(xc),single(zc));
    for ii = 0:Npixels-1
        for jj = 0:Npixels-1
             %temp = fAt3D(2*ii+xstart+off,2*jj+ystart+off,zstart,spline_xsize,spline_ysize,spline_zsize,delta_f,coeff);
             temp = fAt3D(2*ii+xstart+off,2*jj+ystart+off,zstart,spline_xsize,spline_ysize,spline_zsize,delta_f,coeff);
             model = temp*Nphotons+bg;
             output(ii+1,jj+1,kk)=model;
            
        end
    end
    kk
    
end
dipshow(output);
temp1= sum(output,1);
temp2 = sum(temp1,2);
mean(temp2(:))
outputtemp = repmat(output,1,1,10);
outputtemp=single(noise(outputtemp,'poisson',1));

[P_Spline_sim CRLB LL]=GPUmleFit_LM(single(outputtemp(:,:,1:100)),single(coeff*4),50,5,0);

% [P_Spline_sim CRLB LL]=GPUmleFit_LM(single(rand(13,13,10000)),single(coeff*4),50,5,0);
% 
% 
resultfit = [P_Spline_sim(1:100,1:2) P_Spline_sim(1:100,5)];
resultsim = [coordsxy coordsz];
resultfit(1:35,:)-resultsim
% 
% 
% 
% 
% 
% 
% 
% 
% resultfit = [P(:,1:2) P(:,5)];
% resultfit_pixThread = [P_pixThread(:,1:2) P_pixThread(:,5)];
% resultsim = [coordsxy coordsz];
% a1=resultfit(1:100,:)-resultsim;
% a2=resultfit_pixThread(1:100,:)-resultsim;
% 
% 
 %% shift coords

% clear all
% 
% load('matlab161027.mat')

Nfits = 10000;
Nphotons = 5000;
Nphotons =Nphotons*4;
Npixels = 13;
bg = 1;


coordsxy = Npixels/2 -1 +ones([Nfits 2])-0.5;
coordsz = linspace(0,25,Nfits);

% coordsx = Npixels/2 -2 +linspace(0,4,Nfits);
% coordsy = Npixels/2 -1 +ones([Nfits 1])-0.5;
% coordsxy = [coordsx' coordsy];
% coordsz = 8*ones(Nfits,1);

% coordsxy = Npixels/2 -1 +rand([Nfits 2]);
% coordsz = 25*rand(Nfits,1);

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
temp1= sum(output,1);
temp2 = sum(temp1,2);
mean(temp2(:))
outputtemp = repmat(output,1,1,1);
 outputtemp = poissrnd(outputtemp,13,13,100); 
% outputtemp
% outputtemp=single(noise(outputtemp,'poisson',1));

[P_Spline_sim CRLB LL]=GPUmleFit_LM(single(outputtemp(:,:,1:10000)),single(coeff*4),50,5,0);


resultfit = [P_Spline_sim(1:100,1:2) P_Spline_sim(1:100,5)];
resultsim = [coordsxy coordsz];
resultfit(1:100,:)-resultsim








resultfit = [P(:,1:2) P(:,5)];
resultfit_pixThread = [P_pixThread(:,1:2) P_pixThread(:,5)];
resultsim = [coordsxy coordsz];
a1=resultfit(1:100,:)-resultsim;
a2=resultfit_pixThread(1:100,:)-resultsim;


%% with backup j
np_psf = double(psf);
np_psf = np_psf/max(np_psf(:));
% np_psf = np_psf(25-18:25+17,25-18:25+17,:);
np_psf = np_psf(25-12:25+13,25-12:25+13,:);
% np_psf = np_psf(13-6:13+6,13-6:13+6,19-15:19+15);
for i = 1:size(np_psf,3)
    np_psf(:,:,i) = np_psf(:,:,i)/sum(sum(np_psf(:,:,i)));
end
spline = Spline3D(np_psf);
coeff = spline.coeff;

Nfits = 5000;
Nphotons = 5000;
Nphotons =Nphotons*4;
Npixels = 13;
bg = 1;

% coordsx = linspace(Npixels/2-1,Npixels/2+1,Nfits);
% coordsy = linspace(Npixels/2-1,Npixels/2+1,Nfits);
% coordsxy = [coordsx' coordsy'];
% coordsz = linspace(0,25,Nfits);
% coordsxyz = [coordsx' coordsy' coordsz'];
coordsxy = Npixels/2 -3 +6*rand([Nfits 2]);
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
    [delta_f,delta_dxf,delta_ddxf,delta_dyf,delta_ddyf,delta_dzf,delta_ddzf]=computeDelta3Dj(single(xc),single(yc),single(zc));
    for ii = 0:Npixels-1
        for jj = 0:Npixels-1
             temp = fAt3Dj(2*ii+xstart+off,2*jj+ystart+off,zstart,spline_xsize,spline_ysize,spline_zsize,delta_f,coeff);
             model = temp*Nphotons+bg;
             output(ii+1,jj+1,kk)=model;
            
        end
    end
    kk
    
end
temp1= sum(output-bg,1);
temp2 = sum(temp1,2);
mean(temp2(:))
% outputtemp = repmat(output,1,1,1000);
outputtemp=single(noise(output,'poisson',1));

[P_Spline_sim CRLB LL]=GPUmleFit_LM(single(outputtemp(:,:,:)),single(coeff*4),50,5,0);
figure,plot(P_Spline_sim(:,end),'.');
figure,histogram(P_Spline_sim(:,end))
temp =outputtemp(:,:,1:100);
% iw = im2uint16(double(temp));
imwrite(uint16(temp(:,:,1)),'outputtemp.tif');
for ii=2:100
    imwrite(uint16(temp(:,:,ii)),'outputtemp.tif','WriteMode','append');
end

resultfit = [P_Spline_sim(1:100,1:2) P_Spline_sim(1:100,5)];
resultsim = [coordsxy coordsz];
resultfit(1:100,:)-resultsim


%% test stripe
Nfits = 1000;
Nphotons = 1000;
Nphotons =Nphotons*4;
Npixels = 26;
bg = 1;

% coordsx = linspace(Npixels/2-1,Npixels,Nfits);
% coordsy = linspace(Npixels/2-1,Npixels,Nfits);
% coordsz = linspace(0,25,10000);
% coordsxyz = [coordsx' coordsy' coordsz'];
coordsxy = Npixels/2 -3 +6*rand([Nfits 2]);
coordsz = 25*rand(Nfits,1);

output = single(zeros(Npixels,Npixels,Nfits));
spline_xsize = size(coeff,1);
spline_ysize = size(coeff,2);
spline_zsize = size(coeff,3);
off = ((spline_xsize+1)-Npixels);

for kk = 1:Nfits
    xcenter = coordsxy(kk,1);
    ycenter = coordsxy(kk,2);
    zcenter = coordsz(kk);
    xc = (xcenter - 13);
    yc = (ycenter - 13);
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
    [delta_f,delta_dxf,delta_ddxf,delta_dyf,delta_ddyf,delta_dzf,delta_ddzf]=computeDelta3Dj(single(xc),single(yc),single(zc));
    for ii = 0:Npixels-1
        for jj = 0:Npixels-1
             temp = fAt3Dj(ii+xstart+off,jj+ystart+off,zstart,spline_xsize,spline_ysize,spline_zsize,delta_f,coeff);
             model = temp*Nphotons+bg;
             output(ii+1,jj+1,kk)=model;
            
        end
    end
    kk
    
end

output1 = single(zeros(Npixels,Npixels,Npixels));
for kk = 0:Npixels-1
    for ii = 0:Npixels-1
        for jj = 0:Npixels-1
            temp = fAt3Dj(ii,jj,kk,spline_xsize,spline_ysize,spline_zsize,delta_f,coeff);
%             model = temp*Nphotons+bg;
            output1(ii+1,jj+1,kk+1)=temp;
            
        end
    end
end





%%  no interp

%% with backup j
Nfits = 1000;
Nphotons = 1000;
Nphotons =Nphotons*4;
Npixels = 13;
bg = 1;

% coordsx = linspace(Npixels/2-1,Npixels,Nfits);
% coordsy = linspace(Npixels/2-1,Npixels,Nfits);
% coordsz = linspace(0,25,10000);
% coordsxyz = [coordsx' coordsy' coordsz'];
coordsxy = Npixels/2 -3 +6*rand([Nfits 2]);
coordsz = 25*rand(Nfits,1);

output = single(zeros(Npixels,Npixels,Nfits));
spline_xsize = size(coeff1,1);
spline_ysize = size(coeff1,2);
spline_zsize = size(coeff1,3);
off = ((spline_xsize+1)-Npixels);

for kk = 1:Nfits
    xcenter = coordsxy(kk,1);
    ycenter = coordsxy(kk,2);
    zcenter = coordsz(kk);
    xc = (xcenter - 6.5+0.5);
    yc = (ycenter - 6.5+0.5);
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
    [delta_f,delta_dxf,delta_ddxf,delta_dyf,delta_ddyf,delta_dzf,delta_ddzf]=computeDelta3Dj(single(xc),single(yc),single(zc));
    for ii = 0:Npixels-1
        for jj = 0:Npixels-1
             temp = fAt3Dj(ii+xstart+off,jj+ystart+off,zstart,spline_xsize,spline_ysize,spline_zsize,delta_f,coeff1);
             model = temp*Nphotons+bg;
             output(ii+1,jj+1,kk)=model;
            
        end
    end
    kk
    
end
temp1= sum(output-bg,1);
temp2 = sum(temp1,2);
mean(temp2(:))
% outputtemp = repmat(output,1,1,1000);
outputtemp=single(noise(output,'poisson',1));

[P_Spline_sim CRLB LL]=GPUmleFit_LM_noInterp(single(outputtemp(:,:,1:1000)),single(coeff1),50,5,0);
figure,plot(P_Spline_sim(:,end));
figure,histogram(P_Spline_sim(:,end))


resultfit = [P_Spline_sim(:,1:2) P_Spline_sim(:,5)];
resultsim = [coordsxy coordsz];
resultfit-resultsim;
%% back up
Nfits = 100;
Nphotons = 5000;
 Nphotons =Nphotons*4;
Npixels = 13;
bg = 1;

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
    [delta_f,delta_dxf,delta_ddxf,delta_dyf,delta_ddyf,delta_dzf,delta_ddzf]=computeDelta3Dj(single(xc),single(yc),single(zc));
    for ii = 0:Npixels-1
        for jj = 0:Npixels-1
             temp = fAt3Dj(2*ii+xstart+off,2*jj+ystart+off,zstart,spline_xsize,spline_ysize,spline_zsize,delta_f,coeff);
             model = temp*Nphotons+bg;
             output(ii+1,jj+1,kk)=model;
            
        end
    end
    kk
    
end
temp1= sum(output-bg,1);
temp2 = sum(temp1,2);
mean(temp2(:))
% outputtemp = repmat(output,1,1,1000);
% outputtemp=single(noise(output,'poisson',1));
 outputtemp = poissrnd(output,13,13,Nfits); 
[P_Spline_sim CRLB LL]=GPUmleFit_LM(single(outputtemp(:,:,1:10000)),single(coeff*4),50,5,0);
temp =outputtemp(:,:,1:100);
% iw = im2uint16(double(temp));
imwrite(uint16(temp(:,:,1)),'outputtemp.tif');
for ii=2:100
    imwrite(uint16(temp(:,:,ii)),'outputtemp.tif','WriteMode','append');
end

resultfit = [P_Spline_sim(1:100,1:2) P_Spline_sim(1:100,5)];
resultsim = [coordsxy coordsz];
resultfit(1:100,:)-resultsim




%% test Spline3D_v2 and new fitter with increase ROI


np_psf = double(psf);
np_psf = np_psf/max(np_psf(:));
np_psf = np_psf(25-18:25+18,25-18:25+18,:);
% np_psf = np_psf(25-18:25+17,25-18:25+17,:);
% np_psf = np_psf(25-12:25+13,25-12:25+13,:);
% np_psf = np_psf(13-6:13+6,13-6:13+6,19-15:19+15);
for i = 1:size(np_psf,3)
    np_psf(:,:,i) = np_psf(:,:,i)/sum(sum(np_psf(:,:,i)));
end
spline = Spline3D_v2(np_psf);
coeff = spline.coeff;





Nfits = 500;
Nphotons = 5000;
Nphotons =Nphotons*4;
Npixels = 13;
bg = 1;

% coordsx = linspace(Npixels/2-1,Npixels/2+1,Nfits);
% coordsy = linspace(Npixels/2-1,Npixels/2+1,Nfits);
% coordsxy = [coordsx' coordsy'];
% coordsz = linspace(0,25,Nfits);
% coordsxyz = [coordsx' coordsy' coordsz'];
coordsxy = Npixels/2 -3 +6*rand([Nfits 2]);
coordsz = 25*rand(Nfits,1);

output = single(zeros(Npixels,Npixels,Nfits));
spline_xsize = size(coeff,1);
spline_ysize = size(coeff,2);
spline_zsize = size(coeff,3);
% off = ((spline_xsize+1)-2*Npixels)/2;
 off = (spline_xsize-(2*(Npixels)))/2;

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
    [delta_f,delta_dxf,delta_ddxf,delta_dyf,delta_ddyf,delta_dzf,delta_ddzf]=computeDelta3Dj_v2(single(xc),single(yc),single(zc));
    for ii = 0:Npixels-1
        for jj = 0:Npixels-1
             temp = fAt3Dj_v2(2*ii+xstart+off,2*jj+ystart+off,zstart,spline_xsize,spline_ysize,spline_zsize,delta_f,coeff);
             model = temp*Nphotons+bg;
             output(ii+1,jj+1,kk)=model;
            
        end
    end
    kk
    
end
temp1= sum(output-bg,1);
temp2 = sum(temp1,2);
mean(temp2(:))
% outputtemp = repmat(output,1,1,1000);
outputtemp=single(noise(output,'poisson',1));

[P_Spline_sim CRLB LL]=GPUmleFit_LM_v2(single(outputtemp(:,:,:)),single(coeff*4),50,5,0);
figure,plot(P_Spline_sim(:,end),'.');
figure,histogram(P_Spline_sim(:,end))
temp =outputtemp(:,:,1:100);
% iw = im2uint16(double(temp));
imwrite(uint16(temp(:,:,1)),'outputtemp.tif');
for ii=2:100
    imwrite(uint16(temp(:,:,ii)),'outputtemp.tif','WriteMode','append');
end

resultfit = [P_Spline_sim(1:100,1:2) P_Spline_sim(1:100,5)];
resultsim = [coordsxy coordsz];
resultfit(1:100,:)-resultsim





%% test different z plane with Spline3D_v2
 np_psf = double(psf);
np_psf = np_psf/max(np_psf(:));
% np_psf = np_psf(25-12:25+13,25-12:25+13,:);
np_psf = np_psf(25-18:25+17,25-18:25+17,:);
% np_psf = np_psf(13-6:13+6,13-6:13+6,19-15:19+15);
for i = 1:size(np_psf,3)
    np_psf(:,:,i) = np_psf(:,:,i)/sum(sum(np_psf(:,:,i)));
end
spline = Spline3D_v2(np_psf);
coeff = spline.coeff;



 for i = 1:24
Nfits = 1000;
Nphotons = 1000;
Npixels = 13;
bg = 1;

coordsxy = Npixels/2 -1 +rand([Nfits 2]);
% coordsz = 25*rand(Nfits,1);
coordsz = i*ones(Nfits,1);


output = single(zeros(Npixels,Npixels,Nfits));
spline_xsize = size(coeff,1);
spline_ysize = size(coeff,2);
spline_zsize = size(coeff,3);
off = floor(((spline_xsize+1)-2*Npixels)/2);



for kk = 1:Nfits
    xcenter = coordsxy(kk,1);
    ycenter = coordsxy(kk,2);
    zcenter = coordsz(kk);
    
%      zcenter = 1.2;
    %     xc = (xcenter - 6.5+0.5);
    %     yc = (ycenter - 6.5+0.5);
    % for kk = 1:Npixels
    %     xstart = floor(xi);
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
    
   
    [delta_f,delta_dxf,delta_ddxf,delta_dyf,delta_ddyf,delta_dzf,delta_ddzf]=computeDelta3Dj_v2(single(xc),single(yc),single(zc));
    
    for ii = 0:Npixels-1
        for jj = 0:Npixels-1
             temp = fAt3Dj_v2(2*ii+xstart+off,2*jj+ystart+off,zstart,spline_xsize,spline_ysize,spline_zsize,delta_f,coeff);
             model = temp*Nphotons+bg;
             output(ii+1,jj+1,kk)=model;
%             [f,dfx,dfy,dfz] = fAt3Dj_v2(xstart+jj+off,ystart+ii+off,zstart,dataSize(1),dataSize(2),dataSize(3),coeffs,delta_f,delta_dfx,delta_dfy,delta_dfz);
%             %              pd = fAt3D_bs(ii,jj,kk,26,26,26,b3.coeffs);
%             %temp = fAt3D_bs(ii+xstart+off,jj+ystart+off,zstart,spline_xsize,spline_ysize,spline_zsize,delta_f,coeff1);
% %             model = pd*Nphotons+bg;
%             output(ii,jj,kk)=f*Nphotons+bg;
            
        end
    end
%     output(:,:,kk)=output(ii,jj,)
    kk
end

temp1= sum(output-bg,1);
temp2 = sum(temp1,2);
mean(temp2(:))
% 
% outputtemp = repmat(output,1,1,1000);
% outputtemp=single(noise(outputtemp,'poisson',1));
% [P CRLB LL]=GPUmleFit_LM(single(outputtemp(:,:,1:5000)),single(coeff*4),200,5,0);
 output = poissrnd(output,Npixels,Npixels,Nfits); 
% output = single(noise(output,'poisson',1));
[P CRLB LL]=GPUmleFit_LM_v2(single(output),single(coeff),200,5,0);
CRLBx=CRLB(:,1);
CRLBy=CRLB(:,2);
X=P(:,1);
Y=P(:,2);
CRLBz=CRLB(:,5);
Z=P(:,5);
N=P(:,3);
BG=P(:,4);

%report some details
s_x_found(i)=std(X-coordsxy(:,1));
s_y_found(i)=std(Y-coordsxy(:,2));
s_z_found(i)=std(Z-i);

meanCRLBx(i) = mean(sqrt(CRLBx));
meanCRLBy(i) = mean(sqrt(CRLBy));
meanCRLBz(i) = mean(sqrt(CRLBz));

fprintf('The standard deviation of x-position error is %g \n',s_x_found(i))
fprintf('The mean returned CRLB based x-position uncertainty is %g \n',mean(sqrt(CRLBx)))

fprintf('The standard deviation of y-position error is %g \n',s_y_found(i))
fprintf('The mean returned CRLB based y-position uncertainty is %g \n',mean(sqrt(CRLBy)))

fprintf('The standard deviation of z error is %g \n',std(Z-i+1))
fprintf('The mean returned CRLB based z uncertainty is %g \n',mean(sqrt(CRLBz)))
i

 end
 Nz = 24;
s_x_found=s_x_found(1:Nz);
s_y_found=s_y_found(1:Nz);
s_z_found=s_z_found(1:Nz);

meanCRLBx=meanCRLBx(1:Nz);
meanCRLBy=meanCRLBy(1:Nz);
meanCRLBz=meanCRLBz(1:Nz);


s_x_found = s_x_found';
s_y_found = s_y_found';
s_z_found = s_z_found';

meanCRLBx = meanCRLBx';
meanCRLBy = meanCRLBy';
meanCRLBz = meanCRLBz';

%% test different z plane with Spline3D_v2 and odd calibration ROI
 np_psf = double(psf);
np_psf = np_psf/max(np_psf(:));
% np_psf = np_psf(25-12:25+13,25-12:25+13,:);
np_psf = np_psf(25-18:25+17,25-18:25+17,:);
% np_psf = np_psf(13-6:13+6,13-6:13+6,19-15:19+15);
for i = 1:size(np_psf,3)
    np_psf(:,:,i) = np_psf(:,:,i)/sum(sum(np_psf(:,:,i)));
end
spline = Spline3D_v2(np_psf);
coeff = spline.coeff;



 for i = 1:24
Nfits = 1000;
Nphotons = 1000;
Npixels = 13;
bg = 1;

coordsxy = Npixels/2 -1 +rand([Nfits 2]);
% coordsz = 25*rand(Nfits,1);
coordsz = i*ones(Nfits,1);


output = single(zeros(Npixels,Npixels,Nfits));
spline_xsize = size(coeff,1);
spline_ysize = size(coeff,2);
spline_zsize = size(coeff,3);
off = floor(((spline_xsize+1)-2*Npixels)/2);



for kk = 1:Nfits
    xcenter = coordsxy(kk,1);
    ycenter = coordsxy(kk,2);
    zcenter = coordsz(kk);
    
%      zcenter = 1.2;
    %     xc = (xcenter - 6.5+0.5);
    %     yc = (ycenter - 6.5+0.5);
    % for kk = 1:Npixels
    %     xstart = floor(xi);
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
    
   
    [delta_f,delta_dxf,delta_ddxf,delta_dyf,delta_ddyf,delta_dzf,delta_ddzf]=computeDelta3Dj_v2(single(xc),single(yc),single(zc));
    
    for ii = 0:Npixels-1
        for jj = 0:Npixels-1
             temp = fAt3Dj_v2(2*ii+xstart+off,2*jj+ystart+off,zstart,spline_xsize,spline_ysize,spline_zsize,delta_f,coeff);
             model = temp*Nphotons+bg;
             output(ii+1,jj+1,kk)=model;
%             [f,dfx,dfy,dfz] = fAt3Dj_v2(xstart+jj+off,ystart+ii+off,zstart,dataSize(1),dataSize(2),dataSize(3),coeffs,delta_f,delta_dfx,delta_dfy,delta_dfz);
%             %              pd = fAt3D_bs(ii,jj,kk,26,26,26,b3.coeffs);
%             %temp = fAt3D_bs(ii+xstart+off,jj+ystart+off,zstart,spline_xsize,spline_ysize,spline_zsize,delta_f,coeff1);
% %             model = pd*Nphotons+bg;
%             output(ii,jj,kk)=f*Nphotons+bg;
            
        end
    end
%     output(:,:,kk)=output(ii,jj,)
    kk
end

temp1= sum(output-bg,1);
temp2 = sum(temp1,2);
mean(temp2(:))
% 
% outputtemp = repmat(output,1,1,1000);
% outputtemp=single(noise(outputtemp,'poisson',1));
% [P CRLB LL]=GPUmleFit_LM(single(outputtemp(:,:,1:5000)),single(coeff*4),200,5,0);
 output = poissrnd(output,Npixels,Npixels,Nfits); 
% output = single(noise(output,'poisson',1));
[P CRLB LL]=GPUmleFit_LM_v3(single(output),single(coeff1),200,5,0);
CRLBx=CRLB(:,1);
CRLBy=CRLB(:,2);
X=P(:,1);
Y=P(:,2);
CRLBz=CRLB(:,5);
Z=P(:,5);
N=P(:,3);
BG=P(:,4);

%report some details
s_x_found(i)=std(X-coordsxy(:,1));
s_y_found(i)=std(Y-coordsxy(:,2));
s_z_found(i)=std(Z-i);

meanCRLBx(i) = mean(sqrt(CRLBx));
meanCRLBy(i) = mean(sqrt(CRLBy));
meanCRLBz(i) = mean(sqrt(CRLBz));

fprintf('The standard deviation of x-position error is %g \n',s_x_found(i))
fprintf('The mean returned CRLB based x-position uncertainty is %g \n',mean(sqrt(CRLBx)))

fprintf('The standard deviation of y-position error is %g \n',s_y_found(i))
fprintf('The mean returned CRLB based y-position uncertainty is %g \n',mean(sqrt(CRLBy)))

fprintf('The standard deviation of z error is %g \n',std(Z-i+1))
fprintf('The mean returned CRLB based z uncertainty is %g \n',mean(sqrt(CRLBz)))
i

 end
 Nz = 24;
s_x_found=s_x_found(1:Nz);
s_y_found=s_y_found(1:Nz);
s_z_found=s_z_found(1:Nz);

meanCRLBx=meanCRLBx(1:Nz);
meanCRLBy=meanCRLBy(1:Nz);
meanCRLBz=meanCRLBz(1:Nz);


s_x_found = s_x_found';
s_y_found = s_y_found';
s_z_found = s_z_found';

meanCRLBx = meanCRLBx';
meanCRLBy = meanCRLBy';
meanCRLBz = meanCRLBz';

figure,plot(meanCRLBx/2)
hold on
plot(s_x_found)
plot(meanCRLBy/2)

plot(s_y_found)





%% test Spline3D_v2 without double sampling
 np_psf = double(psf);
np_psf = np_psf/max(np_psf(:));
% np_psf = np_psf(25-12:25+13,25-12:25+13,:);
np_psf = np_psf(13-8:13+8,13-8:13+8,:);
% np_psf = np_psf(13-6:13+6,13-6:13+6,19-15:19+15);
for i = 1:size(np_psf,3)
    np_psf(:,:,i) = np_psf(:,:,i)/sum(sum(np_psf(:,:,i)));
end
spline = Spline3D_v2(np_psf);
coeff = spline.coeff;






Nfits = 1000;
Nphotons = 1000;
Nphotons =Nphotons*4;
Npixels = 13;
bg = 1;

% coordsx = linspace(Npixels/2-1,Npixels,Nfits);
% coordsy = linspace(Npixels/2-1,Npixels,Nfits);
% coordsz = linspace(0,25,10000);
% coordsxyz = [coordsx' coordsy' coordsz'];
coordsxy = Npixels/2 -3 +6*rand([Nfits 2]);
coordsz = 25*rand(Nfits,1);

output = single(zeros(Npixels,Npixels,Nfits));
spline_xsize = size(coeff1,1);
spline_ysize = size(coeff1,2);
spline_zsize = size(coeff1,3);
off = ((spline_xsize+1)-Npixels)/2;

for kk = 1:Nfits
    xcenter = coordsxy(kk,1);
    ycenter = coordsxy(kk,2);
    zcenter = coordsz(kk);
     xc = -1*(xcenter - 6.5+0.5);
    yc = -1*(ycenter - 6.5+0.5);
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
    [delta_f,delta_dxf,delta_ddxf,delta_dyf,delta_ddyf,delta_dzf,delta_ddzf]=computeDelta3Dj_v2(single(xc),single(yc),single(zc));
    for ii = 0:Npixels-1
        for jj = 0:Npixels-1
             temp = fAt3Dj_v2(ii+xstart+off,jj+ystart+off,zstart,spline_xsize,spline_ysize,spline_zsize,delta_f,coeff1);
             model = temp*Nphotons+bg;
             output(ii+1,jj+1,kk)=model;
            
        end
    end
    kk
    
end
temp1= sum(output-bg,1);
temp2 = sum(temp1,2);
mean(temp2(:))
% outputtemp = repmat(output,1,1,1000);
outputtemp=single(noise(output,'poisson',1));


% [P] =  kernel_MLEfit_Spline_LM_SMAP_v2_nointerp(single(outputtemp(:,:,1:10)),single(coeff1),13,200);
[P_Spline_sim CRLB LL]=GPUmleFit_LM_noInterp(single(outputtemp(:,:,:)),single(coeff1),200,5,0);

% figure,plot(P_Spline_sim(:,end));
% figure,histogram(P_Spline_sim(:,end))


resultfit = [P_Spline_sim(:,1:2) P_Spline_sim(:,5)];
resultsim = [coordsxy coordsz];
resultfit-resultsim;



%% test different z plane with Spline3D_v2 and odd calibration ROI without double sampling
 np_psf = double(psf);
np_psf = np_psf/max(np_psf(:));
% np_psf = np_psf(25-12:25+13,25-12:25+13,:);
np_psf = np_psf(13-8:13+8,13-8:13+8,:);
% np_psf = np_psf(13-6:13+6,13-6:13+6,19-15:19+15);
for i = 1:size(np_psf,3)
    np_psf(:,:,i) = np_psf(:,:,i)/sum(sum(np_psf(:,:,i)));
end
spline = Spline3D_v2(np_psf);
coeff = spline.coeff;



 for i = 1:24
Nfits = 1000;
Nphotons = 1000;
Npixels = 13;
bg = 1;

coordsxy = Npixels/2 -3 +6*rand([Nfits 2]);
% coordsz = 25*rand(Nfits,1);
coordsz = i*ones(Nfits,1);

output = single(zeros(Npixels,Npixels,Nfits));
spline_xsize = size(coeff,1);
spline_ysize = size(coeff,2);
spline_zsize = size(coeff,3);
off = ((spline_xsize+1)-Npixels)/2;


output = single(zeros(Npixels,Npixels,Nfits));
spline_xsize = size(coeff,1);
spline_ysize = size(coeff,2);
spline_zsize = size(coeff,3);



for kk = 1:Nfits
    xcenter = coordsxy(kk,1);
    ycenter = coordsxy(kk,2);
    zcenter = coordsz(kk);
    
%      zcenter = 1.2;
    %     xc = (xcenter - 6.5+0.5);
    %     yc = (ycenter - 6.5+0.5);
    % for kk = 1:Npixels
    %     xstart = floor(xi);
    xc = -1*(xcenter - 6.5+0.5);
    yc = -1*(ycenter - 6.5+0.5);
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
    
   
    [delta_f,delta_dxf,delta_ddxf,delta_dyf,delta_ddyf,delta_dzf,delta_ddzf]=computeDelta3Dj_v2(single(xc),single(yc),single(zc));
    
    for ii = 0:Npixels-1
        for jj = 0:Npixels-1
             temp = fAt3Dj_v2(ii+xstart+off,jj+ystart+off,zstart,spline_xsize,spline_ysize,spline_zsize,delta_f,coeff);
             model = temp*Nphotons+bg;
             output(ii+1,jj+1,kk)=model;
%             [f,dfx,dfy,dfz] = fAt3Dj_v2(xstart+jj+off,ystart+ii+off,zstart,dataSize(1),dataSize(2),dataSize(3),coeffs,delta_f,delta_dfx,delta_dfy,delta_dfz);
%             %              pd = fAt3D_bs(ii,jj,kk,26,26,26,b3.coeffs);
%             %temp = fAt3D_bs(ii+xstart+off,jj+ystart+off,zstart,spline_xsize,spline_ysize,spline_zsize,delta_f,coeff1);
% %             model = pd*Nphotons+bg;
%             output(ii,jj,kk)=f*Nphotons+bg;
            
        end
    end
%     output(:,:,kk)=output(ii,jj,)
    kk
end

temp1= sum(output-bg,1);
temp2 = sum(temp1,2);
mean(temp2(:))
% 
% outputtemp = repmat(output,1,1,1000);
% outputtemp=single(noise(outputtemp,'poisson',1));
% [P CRLB LL]=GPUmleFit_LM(single(outputtemp(:,:,1:5000)),single(coeff*4),200,5,0);
 output = poissrnd(output,Npixels,Npixels,Nfits); 
% output = single(noise(output,'poisson',1));
[P CRLB LL]=GPUmleFit_LM_noInterp(single(output),single(coeff),200,5,0);
CRLBx=CRLB(:,1);
CRLBy=CRLB(:,2);
X=P(:,1);
Y=P(:,2);
CRLBz=CRLB(:,5);
Z=P(:,5);
N=P(:,3);
BG=P(:,4);

%report some details
s_x_found(i)=std(X-coordsxy(:,1));
s_y_found(i)=std(Y-coordsxy(:,2));
s_z_found(i)=std(Z-i);

meanCRLBx(i) = mean(sqrt(CRLBx));
meanCRLBy(i) = mean(sqrt(CRLBy));
meanCRLBz(i) = mean(sqrt(CRLBz));

fprintf('The standard deviation of x-position error is %g \n',s_x_found(i))
fprintf('The mean returned CRLB based x-position uncertainty is %g \n',mean(sqrt(CRLBx)))

fprintf('The standard deviation of y-position error is %g \n',s_y_found(i))
fprintf('The mean returned CRLB based y-position uncertainty is %g \n',mean(sqrt(CRLBy)))

fprintf('The standard deviation of z error is %g \n',std(Z-i+1))
fprintf('The mean returned CRLB based z uncertainty is %g \n',mean(sqrt(CRLBz)))
i

 end
 Nz = 24;
s_x_found=s_x_found(1:Nz);
s_y_found=s_y_found(1:Nz);
s_z_found=s_z_found(1:Nz);

meanCRLBx=meanCRLBx(1:Nz);
meanCRLBy=meanCRLBy(1:Nz);
meanCRLBz=meanCRLBz(1:Nz);


s_x_found = s_x_found';
s_y_found = s_y_found';
s_z_found = s_z_found';

meanCRLBx = meanCRLBx';
meanCRLBy = meanCRLBy';
meanCRLBz = meanCRLBz';

figure,plot(meanCRLBx)
hold on
plot(s_x_found)
plot(meanCRLBy)

plot(s_y_found)


%% test fitting on the Gaussian Model without double sampling
 np_psf1 = double(psf1);
np_psf1 = np_psf1/max(np_psf1(:));
% np_psf = np_psf(25-12:25+13,25-12:25+13,:);
np_psf1 = np_psf1(13-8:13+8,13-8:13+8,:);
% np_psf = np_psf(13-6:13+6,13-6:13+6,19-15:19+15);
for i = 1:size(np_psf1,3)
    np_psf1(:,:,i) = np_psf1(:,:,i)/sum(sum(np_psf1(:,:,i)));
end
spline1 = Spline3D_v2(np_psf1);
coeff1 = spline1.coeff;

b3_sin=bsarray(double(np_psf1),'lambda',0.4);
coeffs_bs1= b3_sin.coeffs;
coeffs_bs1(end+1,end+1,end+1)=0;
b3_sin.coeffs = coeffs_bs1;


[XX,YY,ZZ]=meshgrid(1:13,1:13,1:26);
[XX1,YY1,ZZ1]=meshgrid(1:0.5:13,1:0.5:13,1:26);
np_psfq = interp3(XX,YY,ZZ,double(np_psf1),XX1,YY1,ZZ1,'spline');

for i = 1:size(np_psfq,3)
    np_psfq(:,:,i) = np_psfq(:,:,i)/sum(sum(np_psfq(:,:,i)));
end
splineq = Spline3D_v2(np_psfq);
coeffq = splineq.coeff;



 np_psf2 = double(psf2);
np_psf2 = np_psf2/max(np_psf2(:));
% np_psf = np_psf(25-12:25+13,25-12:25+13,:);
np_psf2 = np_psf2(25-18:25+18,25-18:25+18,:);
% np_psf = np_psf(13-6:13+6,13-6:13+6,19-15:19+15);
for i = 1:size(np_psf2,3)
    np_psf2(:,:,i) = np_psf2(:,:,i)/sum(sum(np_psf2(:,:,i)));
end
spline2 = Spline3D_v2(np_psf2);
coeff2 = spline2.coeff;








% Nfits=10000         %number of images to fit
% bg=20;               %background fluorescence in photons/pixel/frame
% Nphotons=10000;      %expected photons/frame
% Npixels=13;         %linear size of fit region in pixels.
% PSFsigma=1;         %PSF sigma in pixels
% Ax=0;               %Aberration Terms
% Bx=0;
% Ay=0;
% By=0;
% gamma=.5 ;          %separation of x/y focal planes
% d=0.5;
% z=.6;               %Fixed z position
% fittype=3;
% %   Generate a stack of images
% out = zeros(Npixels,Npixels,Nfits);
% % coordsxyz = zeros(Nfits,3);
% coordsx = linspace(Npixels/2-1,Npixels/2+1,10000);
% coordsy = linspace(Npixels/2-1,Npixels/2+1,10000);
% % coordsxy = Npixels/2 -2 +4*rand([Nfits 2]);
% % coordsx = linspace(Npixels/2-1,Npixels,10000);
% % coordsy = linspace(Npixels/2-1,Npixels,10000);
% coordsz = linspace(-0.6,0.6,10000);
% coordsxyz = [coordsx' coordsy' coordsz'];
% % coordsxyz = [coordsxy coordsz'];
% for i = 1:Nfits
% %     coords=Npixels/2-1+rand([1 2]);
% %     z = 1.2*rand(1,1)-0.6;
% coords = coordsxyz(i,1:2);
% z = coordsxyz(i,3);
% Sx=PSFsigma*sqrt(1+((z-gamma)/d).^2+Ax*((z-gamma)/d).^3+Bx*((z-gamma)/d).^4);
% Sy=PSFsigma*sqrt(1+((z+gamma)/d).^2+Ay*((z+gamma)/d).^3+By*((z+gamma)/d).^4);
% out(:,:,i) = finitegausspsf(Npixels,[Sx Sy],Nphotons,bg,coords);
% % coordsxyz(i,:)=[coords z];
% i
% end
% %   Corrupt with Poisson noise
% % data=single(noise(out,'poisson',1)); %requires DipImage\
% data = poissrnd(out,13,13,Nfits);





Nfits=5917         %number of images to fit
bg=20;               %background fluorescence in photons/pixel/frame
Nphotons=1000;      %expected photons/frame
Npixels=13;         %linear size of fit region in pixels.
PSFsigma=1;         %PSF sigma in pixels
Ax=0;               %Aberration Terms
Bx=0;
Ay=0;
By=0;
gamma=.5 ;          %separation of x/y focal planes
d=0.5;
z=.6;               %Fixed z position
fittype=3;
%   Generate a stack of images
out = zeros(Npixels,Npixels,Nfits);
% coordsxyz = zeros(Nfits,3);
coordsx = linspace(Npixels/2-1,Npixels/2+1,Nfits);
coordsy = linspace(Npixels/2+1,Npixels/2-1,Nfits);
% coordsxy = Npixels/2 -2 +4*rand([Nfits 2]);
% coordsx = linspace(Npixels/2-1,Npixels,10000);
% coordsy = linspace(Npixels/2-1,Npixels,10000);
coordsz = linspace(-0.6,0.6,Nfits);
coordsxyz = [coordsx' coordsy' coordsz'];
% coordsxyz = [coordsxy coordsz'];
for i = 1:Nfits
%     coords=Npixels/2-1+rand([1 2]);
%     z = 1.2*rand(1,1)-0.6;
coords = coordsxyz(i,1:2);
z = coordsxyz(i,3);
Sx=PSFsigma*sqrt(1+((z-gamma)/d).^2+Ax*((z-gamma)/d).^3+Bx*((z-gamma)/d).^4);
Sy=PSFsigma*sqrt(1+((z+gamma)/d).^2+Ay*((z+gamma)/d).^3+By*((z+gamma)/d).^4);
out(:,:,i) = finitegausspsf(Npixels,[Sx Sy],Nphotons,bg,coords);
% coordsxyz(i,:)=[coords z];
i
end
%   Corrupt with Poisson noise
% data=single(noise(out,'poisson',1)); %requires DipImage\
data = poissrnd(out,13,13,Nfits);




[P_noInterp CRLB LL]=GPUmleFit_LM_noInterp(single(data),single(coeff1),200,5,0);
[P_noInterp CRLB LL]=GPUmleFit_LM_noInterp_cuda8(single(data),single(coeff1),200,5,0);
figure,plot(P_noInterp(:,end),'.')
figure,plot(P_noInterp(:,1),'.')
figure,plot(P_noInterp(:,2),'.')

% temp = P_noInterp(:,1);
% P_noInterp(:,1)=P_noInterp(:,2);
% P_noInterp(:,2)=temp;




[P_v2 CRLB LL]=GPUmleFit_LM_v2(single(data),single(coeff2*4),200,5,0);
figure,plot(P_v2(:,end),'.')
figure,plot(P_v2(:,1),'.')
figure,plot(P_v2(:,2),'.')



[P_v2q CRLB LL]=GPUmleFit_LM_v2(single(data),single(coeffq*4),200,5,0);
figure,plot(P_v2q(:,end),'.')
figure,plot(P_v2q(:,1),'.')
figure,plot(P_v2q(:,2),'.')


 
% [P_pix_bs CRLB LL]=GPUmleFit_LM_pixThread_bSpline(single(data),single(coeffs_bs1),200,5,0);
% figure,plot(P_pix_bs(:,end),'.')
% figure,plot(P_pix_bs(:,1),'.')
% figure,plot(P_pix_bs(:,2),'.')
% 
% temp = P_pix_bs(:,1);
% P_pix_bs(:,1)=P_pix_bs(:,2);
% P_pix_bs(:,2)=temp;

[P_pix_bs2 CRLB LL]=GPUmleFit_LM_pixThread_bSpline_v2(single(data),single(coeffs_bs1),200,5,0);
figure,plot(P_pix_bs2(:,end),'.')
figure,plot(P_pix_bs2(:,1),'.')
figure,plot(P_pix_bs2(:,2),'.')


% [P_bs CRLB LL]=GPUmleFit_LM_bSpline(single(data),single(coeffs_bs1),200,5,0);
% figure,plot(P_bs(:,end),'.')
% figure,plot(P_bs(:,1),'.')
% figure,plot(P_bs(:,2),'.')

[P_bs2 CRLB LL]=GPUmleFit_LM_bSpline_v2(single(data),single(coeffs_bs1),200,5,0);
figure,plot(P_bs2(:,end),'.')
figure,plot(P_bs2(:,1),'.')
figure,plot(P_bs2(:,2),'.')




% temp = P_CS1-P_CS2;
% figure,plot(temp(:,end),'.')
% figure,plot(temp(:,1),'.')
% figure,plot(temp(:,2),'.')




% [P_M,update] =  kernel_MLEfit_Spline_LM_v8_single(data(:,:,1:10000),b3_sin,13,200);
% figure,plot(P_M(:,5),'.')
% figure,plot(P_M(:,1),'.')
% figure,plot(P_M(:,2),'.')

temp_bs_bs = P_bs2-P_pix_bs2;
figure,plot(temp_bs_bs(:,end),'.')
figure,plot(temp_bs_bs(:,1),'.')
figure,plot(temp_bs_bs(:,2),'.')


temp_bs_noInterp=P_pix_bs2-P_noInterp;
figure,plot(temp_bs_noInterp(:,end),'.')
figure,plot(temp_bs_noInterp(:,1),'.')
figure,plot(temp_bs_noInterp(:,2),'.')


tempCS_noInterp_v2=P_noInterp-P_v2;
figure,plot(tempCS_noInterp_v2(:,end),'.')
figure,plot(tempCS_noInterp_v2(:,1),'.')
figure,plot(tempCS_noInterp_v2(:,2),'.')


tempCS_noInterp_v2q=P_noInterp-P_v2q;
figure,plot(tempCS_noInterp_v2q(:,end),'.')
figure,plot(tempCS_noInterp_v2q(:,1),'.')
figure,plot(tempCS_noInterp_v2q(:,2),'.')















