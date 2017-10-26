%% test different z plane with Spline3D_v2 and odd calibration ROI without double sampling
%  np_psf = double(psf);
% np_psf = np_psf/max(np_psf(:));
% % np_psf = np_psf(25-12:25+13,25-12:25+13,:);
% np_psf = np_psf(13-8:13+8,13-8:13+8,:);
% % np_psf = np_psf(13-6:13+6,13-6:13+6,19-15:19+15);
% for i = 1:size(np_psf,3)
%     np_psf(:,:,i) = np_psf(:,:,i)/sum(sum(np_psf(:,:,i)));
% end
% spline = Spline3D_v2(np_psf);
% coeff = spline.coeff;

load('10_zstack_1_3Dcal.mat');
coeff = SXY.splinefit.cspline.coeff;
i = 1;

Nfits = 1000;
Nphotons =5000;
Npixels = 13;
bg = 1;

coordsxy = Npixels/2 -1 +2*rand([Nfits 2]);
% coordsz = 25*rand(Nfits,1);
coordsz = n*ones(Nfits,1);

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
    kk;
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
[P CRLB LL]=GPUmleFit_LM_faster(single(output),5,100,single(coeff),0,0);
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


[P1 CRLB1 LL1]=GPUmleFit_LM_New12(single(output),5,100,single(coeff),0,0);
CRLBx1=CRLB1(:,1);
CRLBy1=CRLB1(:,2);
X1=P1(:,1);
Y1=P1(:,2);
CRLBz1=CRLB1(:,5);
Z1=P1(:,5);
N1=P1(:,3);
BG1=P1(:,4);

%report some details
s_x_found1(i)=std(X1-coordsxy(:,1));
s_y_found1(i)=std(Y1-coordsxy(:,2));
s_z_found1(i)=std(Z1-i);

meanCRLBx1(i) = mean(sqrt(CRLBx1));
meanCRLBy1(i) = mean(sqrt(CRLBy1));
meanCRLBz1(i) = mean(sqrt(CRLBz1));

fprintf('The standard deviation of x-position error is %g \n',s_x_found1(i))
fprintf('The mean returned CRLB based x-position uncertainty is %g \n',mean(sqrt(CRLBx1)))

fprintf('The standard deviation of y-position error is %g \n',s_y_found1(i))
fprintf('The mean returned CRLB based y-position uncertainty is %g \n',mean(sqrt(CRLBy1)))

fprintf('The standard deviation of z error is %g \n',std(Z1-i+1))
fprintf('The mean returned CRLB based z uncertainty is %g \n',mean(sqrt(CRLBz1)))

[P2 CRLB2 LL2]=GPUmleFit_LM_Newfit(single(output),5,100,single(coeff),0,0);
CRLBx2=CRLB2(:,1);
CRLBy2=CRLB2(:,2);
X2=P2(:,1);
Y2=P2(:,2);
CRLBz2=CRLB2(:,5);
Z2=P2(:,5);
N2=P2(:,3);
BG2=P2(:,4);

%report some details
s_x_found2(i)=std(X2-coordsxy(:,1));
s_y_found2(i)=std(Y2-coordsxy(:,2));
s_z_found2(i)=std(Z2-i);

meanCRLBx2(i) = mean(sqrt(CRLBx2));
meanCRLBy2(i) = mean(sqrt(CRLBy2));
meanCRLBz2(i) = mean(sqrt(CRLBz2));

fprintf('The standard deviation of x-position error is %g \n',s_x_found2(i))
fprintf('The mean returned CRLB based x-position uncertainty is %g \n',mean(sqrt(CRLBx2)))

fprintf('The standard deviation of y-position error is %g \n',s_y_found2(i))
fprintf('The mean returned CRLB based y-position uncertainty is %g \n',mean(sqrt(CRLBy2)))

fprintf('The standard deviation of z error is %g \n',std(Z2-i+1))
fprintf('The mean returned CRLB based z uncertainty is %g \n',mean(sqrt(CRLBz2)))


%  Nz = 130;
% s_x_found=s_x_found(1:Nz);
% s_y_found=s_y_found(1:Nz);
% s_z_found=s_z_found(1:Nz);
% 
% meanCRLBx=meanCRLBx(1:Nz);
% meanCRLBy=meanCRLBy(1:Nz);
% meanCRLBz=meanCRLBz(1:Nz);


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

plot(s_z_found)
plot(meanCRLBz)


s_x_found1 = s_x_found1';
s_y_found1 = s_y_found1';
s_z_found1 = s_z_found1';

meanCRLBx1 = meanCRLBx1';
meanCRLBy1 = meanCRLBy1';
meanCRLBz1 = meanCRLBz1';

figure,plot(meanCRLBx1)
hold on
plot(s_x_found1)
plot(meanCRLBy1)

plot(s_y_found1)

plot(s_z_found1)
plot(meanCRLBz1)





s_x_found2 = s_x_found2';
s_y_found2 = s_y_found2';
s_z_found2 = s_z_found2';

meanCRLBx2 = meanCRLBx2';
meanCRLBy2 = meanCRLBy2';
meanCRLBz2 = meanCRLBz2';

figure,plot(meanCRLBx2)
hold on
plot(s_x_found2)
plot(meanCRLBy2)

plot(s_y_found2)

plot(s_z_found2)
plot(meanCRLBz2)




