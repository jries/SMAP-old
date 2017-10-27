
%parameters affine transform
% theta = 5*pi/180;%rotation
% sx = 1;%scalex
% sy = 1.1;%scaley
% tx = 1;%shfitx
% ty = 1;%shfity



theta = 0*pi/180;%rotation
sx = 1;%scalex
sy = 1;%scaley
tx = 0.3;%shfitx
ty = 0.2;%shfity

%define affine transformation matrix
tformR = [cos(theta) -sin(theta) 0;sin(theta) cos(theta) 0;0 0 1];
tformS = [sx 0 0;0 sy 0;0 0 1];
tformT = [1 0 tx; 0 1 ty; 0 0 1];


tform = tformR*tformS*tformT;


cal1=load('bead_3dcal_old.mat');
coeff1 = cal1.cspline.coeff;
dz1=cal1.cspline.dz;
z01=cal1.cspline.z0; 


cal2=load('bead_3dcal_old.mat');
coeff2 = cal2.cspline.coeff;
dz2=cal2.cspline.dz;
z02=cal2.cspline.z0; 





Nfits = 1000;
Nphotons1 =5000;
Nphotons2 =3000;
Npixels = 13;
bg1 = 15;
bg2 = 10;

coordsxy1 = Npixels/2 -1 +2*rand([Nfits 2]);
ztruth = 300;
coordsz1 = ztruth/dz1+z01*ones(Nfits,1);
coordsz2 = coordsz1;
% coordsz2 = ztruth/dz2+z02*ones(Nfits,1);

coordsxy2 = zeros(Nfits,2);

% tCoordsxy = tform*[coordsxy';1];
% tCoordsxy = tCoordsxy(1:2);

for i = 1:Nfits
    temp = tform*[coordsxy1(i,:)';1];
    coordsxy2(i,:) = temp(1:2);
end









output1 = single(zeros(Npixels,Npixels,Nfits));
output2 = single(zeros(Npixels,Npixels,Nfits));
spline_xsize1 = size(coeff1,1);
spline_ysize1 = size(coeff1,2);
spline_zsize1 = size(coeff1,3);
off1 = ((spline_xsize1+1)-Npixels)/2;

spline_xsize2 = size(coeff2,1);
spline_ysize2 = size(coeff2,2);
spline_zsize2 = size(coeff2,3);
off2 = ((spline_xsize2+1)-Npixels)/2;



% output = single(zeros(Npixels,Npixels,Nfits));



for kk = 1:Nfits
    kk
    xcenter1 = coordsxy1(kk,1);
    ycenter1 = coordsxy1(kk,2);
    
    xcenter2 = coordsxy2(kk,1);
    ycenter2 = coordsxy2(kk,2);
    
    
    zcenter1 = coordsz1(kk);
    zcenter2 = coordsz2(kk);
%      zcenter = 1.2;
    %     xc = (xcenter - 6.5+0.5);
    %     yc = (ycenter - 6.5+0.5);
    % for kk = 1:Npixels
    %     xstart = floor(xi);
    xc1 = -1*(xcenter1 - Npixels/2+0.5);
    yc1 = -1*(ycenter1 - Npixels/2+0.5);
    
    xc2 = -1*(xcenter2 - Npixels/2+0.5);
    yc2 = -1*(ycenter2 - Npixels/2+0.5);
    
    zc1 = zcenter1 -floor(zcenter1);
    zstart1 = floor(zcenter1);
    
     zc2 = zcenter2 -floor(zcenter2);
    zstart2 = floor(zcenter2);
    
    xstart1 = floor(xc1);
    xc1 = xc1 - xstart1;
    
    ystart1 = floor(yc1);
    yc1 = yc1 - ystart1;
    
    xstart2 = floor(xc2);
    xc2 = xc2 - xstart2;
    
    ystart2 = floor(yc2);
    yc2 = yc2 - ystart2;
    
    
    
    
  
    
   
    [delta_f1,delta_dxf1,delta_ddxf1,delta_dyf1,delta_ddyf1,delta_dzf1,delta_ddzf1]=computeDelta3Dj_v2(single(xc1),single(yc1),single(zc1));
    
    [delta_f2,delta_dxf2,delta_ddxf2,delta_dyf2,delta_ddyf2,delta_dzf2,delta_ddzf2]=computeDelta3Dj_v2(single(xc2),single(yc2),single(zc2));
    
    for ii = 0:Npixels-1
        for jj = 0:Npixels-1
             temp1 = fAt3Dj_v2(ii+xstart1+off1,jj+ystart1+off1,zstart1,spline_xsize1,spline_ysize1,spline_zsize1,delta_f1,coeff1);
             temp2 = fAt3Dj_v2(ii+xstart2+off2,jj+ystart2+off2,zstart2,spline_xsize2,spline_ysize2,spline_zsize2,delta_f2,coeff2);
             model1 = temp1*Nphotons1+bg1;
              model2 = temp2*Nphotons2+bg2;
             output1(ii+1,jj+1,kk)=model1;
             
             output2(ii+1,jj+1,kk)=model2;
%             [f,dfx,dfy,dfz] = fAt3Dj_v2(xstart+jj+off,ystart+ii+off,zstart,dataSize(1),dataSize(2),dataSize(3),coeffs,delta_f,delta_dfx,delta_dfy,delta_dfz);
%             %              pd = fAt3D_bs(ii,jj,kk,26,26,26,b3.coeffs);
%             %temp = fAt3D_bs(ii+xstart+off,jj+ystart+off,zstart,spline_xsize,spline_ysize,spline_zsize,delta_f,coeff1);
% %             model = pd*Nphotons+bg;
%             output(ii,jj,kk)=f*Nphotons+bg;
            
        end
    end
%     output(:,:,kk)=output(ii,jj,)
    kk;
end

% temp1= sum(output-bg,1);
% temp2 = sum(temp1,2);
% mean(temp2(:))
% 
% outputtemp = repmat(output,1,1,1000);
% outputtemp=single(noise(outputtemp,'poisson',1));
% [P CRLB LL]=GPUmleFit_LM(single(outputtemp(:,:,1:5000)),single(coeff*4),200,5,0);
 output1 = poissrnd(output1,Npixels,Npixels,Nfits); 
output2 = poissrnd(output2,Npixels,Npixels,Nfits); 

temp = [output1 output2];
dipshow(temp)
%% fit molecules
coords1chip = zeros(Nfits,2);
coords2chip = zeros(Nfits,2);

iterations=50;
% dxy=tform*coords1chip-coords2chip

[P0,update0, error0, model_save0] =  kernel_MLEfit_Spline_LM_multichannel(output1(:,:,1:10),output2(:,:,1:10),coords1chip, coords2chip, tform, coeff1,coeff2,Npixels,iterations);






%% fit molecules_finalized


noChannels = 2;
iterations=50;
dT = zeros(5,noChannels,Nfits);
dxy=coordsxy1-coordsxy2;
temp = reshape(dxy',[2 1,Nfits]);
dT(1:2,2,:)=temp*-1;
shared = [1;1;1;0;0];
sharedA = repmat(shared,[1 Nfits]);
d_data(:,:,:,1) = output1;
d_data(:,:,:,2) = output2;

coeff(:,:,:,:,1)=coeff1;
coeff(:,:,:,:,2)=coeff2;


% [P,CRLB, LL,update, error] =  kernel_MLEfit_Spline_LM_multichannel_finalized(d_data(:,:,1:10,:),coeff, shared,dT,1);
% 
% 
% [P1,CRLB1, LL1] =  GPUmleFit_LM_MultiChannel(single(d_data(:,:,:,:)),int32(sharedA),50,single(coeff),single(dT));

d_data = single(d_data);
sharedA = int32(sharedA);
coeff = single(coeff);
dT = single(dT);
tic
[P1,CRLB1, LL1] =  GPUmleFit_LM_MultiChannel(d_data,sharedA,50,coeff,dT);

tspline=toc;
disp(['cspline: ' num2str(Nfits/tspline) ' fits/s']);





