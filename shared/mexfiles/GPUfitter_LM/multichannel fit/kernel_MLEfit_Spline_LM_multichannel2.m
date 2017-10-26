function [P,update, error, model_save] =  kernel_MLEfit_Spline_LM_multichannel2(d_data1,d_data2,coords1chip, coords2chip, tform, coeff1,coeff2,sz,iterations)
% update computeDelta3Dj and fAt3Dj. Add kernel_derivativeSpline
% make it capable of handling coeff with random size using Spline3D_v2
%work with odd number of ROI for calibration
%try to solve the negative value problem
pi = single(3.141592);
spline_xsize = size(coeff1, 1);
spline_ysize = size(coeff1, 2);
spline_zsize = size(coeff1, 3);

spline_xsize2 = size(coeff2, 1);
spline_ysize2 = size(coeff2, 2);
spline_zsize2 = size(coeff2,3);


PSFSigma = single(1.5);
xc1 = single(0);
yc1 = single(0);
zc1 = single(0);

xc2 = single(0);
yc2 = single(0);
zc2 = single(0);


NV = 7; %number of parameters to optimize

dudt = single(zeros(sz,sz,NV));

newTheta1 = single(zeros(5,1));
newTheta2 = single(zeros(5,1));

% oldTheta = newTheta;

M = zeros(NV*NV,1);
Minv = zeros(NV*NV,1);

Nfits = size(d_data1,3);

error = zeros(Nfits, 1);

RUNNING = 0;
CONVERGED = 1;
CHOLERRER = 2;
BADPEAK = 3;
tolerance=1e-6;

num_of_iterations = 0;
model_save = [];
sigma_x = zeros(Nfits, 1);
sigma_y = zeros(Nfits, 1);
sigma_to_z = 0.55;

weightPara = [1,1;1,1;1,1;1,0;0,1;1,0;0,1];
for tx = 1:Nfits
    tx;
    %newpeak
    [newTheta(1),newTheta(2)]=kernel_CenterofMass2D(sz,single(d_data1(sz*sz*(tx-1)+1:sz*sz*(tx))));
    temp = tform*[newTheta(1)+coords1chip(tx,1);newTheta(2)+coords1chip(tx,2);1];
    newThetaT(1) = temp(1)-coords2chip(tx,1);
    newThetaT(2) = temp(2)-coords2chip(tx,2);
    
    [Nmax1,newTheta(6)] = kernel_GaussFMaxMin2D(sz,PSFSigma,single(d_data1(sz*sz*(tx-1)+1:sz*sz*(tx))));  
    [Nmax2,newTheta(7)] = kernel_GaussFMaxMin2D(sz,PSFSigma,single(d_data2(sz*sz*(tx-1)+1:sz*sz*(tx))));
%     newTheta(7) = newTheta(6);
    newTheta(4)=max(0,(Nmax1-newTheta(6))*2*pi*PSFSigma*PSFSigma);
    newTheta(5)=max(0,(Nmax2-newTheta(7))*2*pi*PSFSigma*PSFSigma);
%     newTheta(5)= newTheta(4);
    newTheta(3) = spline_zsize/2;
%     newTheta2(4) =  newTheta1(4);
%     newTheta(5) = in_pos;
%     newTheta(6) = 0;
    
    newLambda = (0.01);
    newUpdate = 100000*ones(NV,1);
    oldUpdate = 100000*ones(NV,1);
    
    maxJump = [1, 1,spline_zsize/5,100,100,100,100];
    maxJump(4) = max(maxJump(4),newTheta(4));
    maxJump(5) = max(maxJump(5),newTheta(5));
    
    maxJump(6) = max(maxJump(6),newTheta(6));
    maxJump(7) = max(maxJump(7),newTheta(7));
    
%     maxJump(6) = max(maxJump(6),newTheta(6));
%     maxJump(7) = max(maxJump(7),newTheta(7));
    
    
    
    
    %updateFitValues3D
%     xc = single(-2*(newTheta(1) - sz/2+0.5));
%     yc = single(-2*(newTheta(2) - sz/2+0.5));
    xc1 = single(-1 * (newTheta(1) - sz / 2 + 0.5));
    yc1 = single(-1 * (newTheta(2) - sz / 2 + 0.5));
    
    
     xc2 = single(-1 * (newThetaT(1) - sz / 2 + 0.5));
    yc2 = single(-1 * (newThetaT(2) - sz / 2 + 0.5));
    
    
    
    
    zc = single(newTheta(3)-floor(newTheta(3)));
    off = floor(((spline_xsize+1)-sz)/2);
    
    xstart1 = floor(xc1);
    xc1 = xc1-floor(xc1);
    
    ystart1 = floor(yc1);
    yc1 = yc1-floor(yc1);
    
    
    
    xstart2 = floor(xc2);
    xc2 = xc2-floor(xc2);
    
   ystart2 = floor(yc2);
    yc2 = yc2-floor(yc2);
    
    
    
    newErr1 = 0;
    newErr2 = 0;
    newErr = 0;
    jacobian = single(zeros(NV, 1));
    hessian = single(zeros(NV, NV));
    zstart = single(floor(newTheta(3)));
    [delta_f1,delta_dxf1,delta_ddxf1,delta_dyf1,delta_ddyf1,delta_dzf1,delta_ddzf1] = computeDelta3Dj_v2(single(xc1),single(yc1),single(zc));
    
     [delta_f2,delta_dxf2,delta_ddxf2,delta_dyf2,delta_ddyf2,delta_dzf2,delta_ddzf2] = computeDelta3Dj_v2(single(xc2),single(yc2),single(zc));
    
     newTheta1 = [newTheta(1),newTheta(2),newTheta(3),newTheta(4),newTheta(6)];
     newTheta2 = [newThetaT(1),newThetaT(2),newTheta(3),newTheta(5),newTheta(7)];
    for ii = 0:sz-1
        for jj = 0:sz-1
            
            data1 = single(d_data1(sz*sz*(tx-1)+sz*jj+ii+1));
            data2 = single(d_data2(sz*sz*(tx-1)+sz*jj+ii+1));
            
            
            [newDudt1, model1] =  kernel_DerivativeSpline_v21(ii+xstart1+off,jj+ystart1+off,zstart,spline_xsize,spline_ysize,spline_zsize,delta_f1,delta_dxf1,delta_dyf1,delta_dzf1,coeff1,newTheta1,1);
            
            [newDudt2, model2] =  kernel_DerivativeSpline_v21(ii+xstart2+off,jj+ystart2+off,zstart,spline_xsize,spline_ysize,spline_zsize,delta_f2,delta_dxf2,delta_dyf2,delta_dzf2,coeff2,newTheta2,2);
            
            
            
            
%             if data>0
%                 newErr = newErr +2*((model-data)-data*log(model/data));
%             elseif data ==0
%                 newErr = newErr + 2*model;
%             end
            
            if data1>0
                newErr = newErr +2*((model1-data1)-data1*log(model1/data1));
            else
                newErr = newErr + 2*model1;
                data1 = 0;
            end
            
            if data2>0
                newErr = newErr +2*((model2-data2)-data2*log(model2/data2));
            else
                newErr = newErr + 2*model2;
                data2 = 0;
            end
            
            if ~isreal(newErr1)&&~isreal(newErr1)
                ii
                jj
            end
%             newErr1t(ii+1,jj+1)=newErr1;
%             newErr2t(ii+1,jj+1)=newErr2;
%             
%             
%             newErr = newErr+newErr1+newErr2;
            
            
            
            t11 = single(1 - data1 / model1);
            t12 = single(1 - data2 / model2);
            
             for l = 1:NV
                jacobian(l) = jacobian(l)+t11*newDudt1(l)+t12*newDudt2(l);
%                 jacobian2(l) = jacobian2(l)+t12*newDudt2(l);
             end
            
             
            t21 = data1/model1^2;
            t22 = data2/model2^2;
            for l = 0:NV-1
                for m =l:NV-1
                    hessian(l*NV+m+1) = hessian(l*NV+m+1)+t21*newDudt1(l+1)*newDudt1(m+1)+t22*newDudt2(l+1)*newDudt2(m+1);
                    hessian(m*NV+l+1) = hessian(l*NV+m+1);
                end
            end
            
        end
    end
    %         dipshow(newDudt1);
    %     newErr = err;
    oldErr = 1e13;
    %addPeak
    
    %copyFitData
    oldLambda = newLambda;
    oldTheta = newTheta;
    info = 0;
    for kk =1:iterations
        
        
        if abs((newErr-oldErr)/newErr)<tolerance
            %                 newStatus = CONVERGED;
            disp('converged')
            error(tx) = newErr;
            num_of_iterations = kk;
            break;
        else
            if newErr>1.5*oldErr
                newLambda = oldLambda;
                newTheta = oldTheta;
                newErr = oldErr;
                newUpdate = oldUpdate;
%                 if newLambda<0.1
%                     newLambda = 0.1;
%                 end
%                 mult = (1 + newLambda*10)/(1 + newLambda);
%                 maxJump = maxJump*2;
                mult = max( (1 + newLambda*10)/(1 + newLambda),1.3);
                newLambda = 10*newLambda;
                %addpeak
%                 maxJump = maxJump*2;
                
            elseif newErr<oldErr &&info == 0
                
%                     newLambda = 0.1*newLambda;
                    mult = 1 + newLambda;
                    newLambda = 0.1*newLambda;
               
            end
            
            
            for i = 0:NV-1
%                 mult=1
                hessian(i*NV+i+1) = hessian(i*NV+i+1)*mult;
            end
            
            
            
            
            [L, U ,info] = kernel_cholesky(hessian,NV);
            if info ==0
                
                 kk;
                oldLambda = newLambda;
                oldTheta = newTheta;
                 oldUpdate = newUpdate;
                oldErr = newErr;
                
                
                newUpdate = kernel_luEvaluate(L,U,jacobian,NV);
               
                
                %copyFitData
               
               
                
                
                %updatePeakParameters
                
                for ll =1:NV
                    ll;
%                     newUpdate(ll) = oldUpdate(ll)/(1+abs(oldUpdate(ll)/maxJump(ll)));
 %                    temp = min(max(newUpdate(ll),-maxJump(ll)),maxJump(ll));
                    if newUpdate(ll)/oldUpdate(ll)<-0.5
%                         newUpdate(ll) = temp*0.5;
                         maxJump(ll) = maxJump(ll)*0.5;
                         newUpdate(ll) = newUpdate(ll)/(1+abs(newUpdate(ll)/maxJump(ll)));
                    else
                        newUpdate(ll) =  newUpdate(ll)/(1+abs(newUpdate(ll)/maxJump(ll)));
                    end
                    newTheta(ll) = newTheta(ll)-newUpdate(ll);
                    update(kk,ll) = newUpdate(ll);
                end
                 
%                 update(kk,NV+1)= oldErr;
                
                newTheta(1) = max(newTheta(1),(sz-1)/2-sz/4.0);
                newTheta(1) = min(newTheta(1),(sz-1)/2+sz/4.0);
                newTheta(2) = max(newTheta(2),(sz-1)/2-sz/4.0);
                newTheta(2) = min(newTheta(2),(sz-1)/2+sz/4.0);
                newTheta(4) = max(newTheta(4), 1);
                newTheta(5) = max(newTheta(5), 1);
%                 newTheta(6) = max(newTheta(6), 0.01);
%                 newTheta(7) = max(newTheta(7), 0.01);
                newTheta(3) = max(newTheta(3), 0);
                newTheta(3) = min(newTheta(3), spline_zsize);
                
%                 newTheta(6) = max(newTheta(6), 0);
%                 newTheta(6) = min(newTheta(6), 2 * pi);
                
                %updateFitValues3D
                %                 xc = single(-2*((newTheta(1) - sz/2)+0.5));
                %                 yc = single(-2*((newTheta(2) - sz/2)+0.5));
                %                 zc = single(newTheta(5)-floor(newTheta(5)));
                
                
                temp = tform*[newTheta(1)+coords1chip(tx,1);newTheta(2)+coords1chip(tx,2);1];
                newThetaT(1) = temp(1)-coords2chip(tx,1);
                newThetaT(2) = temp(2)-coords2chip(tx,2);
                
                xc1 = single(-1 * (newTheta(1) - sz / 2 + 0.5));
                yc1 = single(-1 * (newTheta(2) - sz / 2 + 0.5));
                
                
                xc2 = single(-1 * (newThetaT(1) - sz / 2 + 0.5));
                yc2 = single(-1 * (newThetaT(2) - sz / 2 + 0.5));
                
                
                
                
                zc = single(newTheta(3)-floor(newTheta(3)));
                off = floor(((spline_xsize+1)-sz)/2);
                
                xstart1 = floor(xc1);
                xc1 = xc1-floor(xc1);
                
                ystart1 = floor(yc1);
                yc1 = yc1-floor(yc1);
                
                
                
                xstart2 = floor(xc2);
                xc2 = xc2-floor(xc2);
                
                ystart2 = floor(yc2);
                yc2 = yc2-floor(yc2);
                
                
                
                newErr1 = 0;
                newErr2 = 0;
                newErr = 0;
                jacobian = single(zeros(NV, 1));
                hessian = single(zeros(NV, NV));
                
                
                %                 xc = single(-1*(newTheta(1) - sz/2+0.5));
                %                 yc = single(-1*(newTheta(2) - sz/2+0.5));
                %                 zc = single(newTheta(5)-floor(newTheta(5)));
                %
                %                 xstart = floor(xc);
                %                 xc = xc-floor(xc);
                %
                %                 ystart = floor(yc);
                %                 yc = yc-floor(yc);
                %
                %                 newErr = 0;
                %                 jacobian = single(zeros(NV,1));
                %                 hessian = single(zeros(NV,NV));
                zstart = single(floor(newTheta(3)));
                %                 [delta_f,delta_dxf,delta_ddxf,delta_dyf,delta_ddyf,delta_dzf,delta_ddzf]=computeDelta3Dj_v2(single(xc),single(yc),single(zc));
                [delta_f1,delta_dxf1,delta_ddxf1,delta_dyf1,delta_ddyf1,delta_dzf1,delta_ddzf1] = computeDelta3Dj_v2(single(xc1),single(yc1),single(zc));
                
                [delta_f2,delta_dxf2,delta_ddxf2,delta_dyf2,delta_ddyf2,delta_dzf2,delta_ddzf2] = computeDelta3Dj_v2(single(xc2),single(yc2),single(zc));
                newTheta1 = [newTheta(1),newTheta(2),newTheta(3),newTheta(4),newTheta(5)];
                newTheta2 = [newThetaT(1),newThetaT(2),newTheta(3),newTheta(4),newTheta(5)];
                for ii = 0:sz-1
                    for jj = 0:sz-1
                        data1 = single(d_data1(sz*sz*(tx-1)+sz*jj+ii+1));
                        data2 = single(d_data2(sz*sz*(tx-1)+sz*jj+ii+1));
                        
                        
                        [newDudt1, model1] =  kernel_DerivativeSpline_v21(ii+xstart1+off,jj+ystart1+off,zstart,spline_xsize,spline_ysize,spline_zsize,delta_f1,delta_dxf1,delta_dyf1,delta_dzf1,coeff1,newTheta1,1);
                        
                        [newDudt2, model2] =  kernel_DerivativeSpline_v21(ii+xstart2+off,jj+ystart2+off,zstart,spline_xsize,spline_ysize,spline_zsize,delta_f2,delta_dxf2,delta_dyf2,delta_dzf2,coeff2,newTheta2,2);
                        
                        
                        
                        
                        %             if data>0
                        %                 newErr = newErr +2*((model-data)-data*log(model/data));
                        %             elseif data ==0
                        %                 newErr = newErr + 2*model;
                        %             end
                        
                        if data1>0
                            newErr = newErr +2*((model1-data1)-data1*log(model1/data1));
                        else
                            newErr = newErr + 2*model1;
                            data1 = 0;
                        end
                        
                        if data2>0
                            newErr = newErr +2*((model2-data2)-data2*log(model2/data2));
                        else
                            newErr = newErr + 2*model2;
                            data2 = 0;
                        end
                        
%                         newErr = newErr+newErr1+newErr2;
                        
                        
                        
                        t11 = single(1 - data1 / model1);
                        t12 = single(1 - data2 / model2);
                        
                        for l = 1:NV
                            jacobian(l) = jacobian(l)+t11*newDudt1(l)+t12*newDudt2(l);
                            %                 jacobian2(l) = jacobian2(l)+t12*newDudt2(l);
                        end
                        
                        
                        t21 = data1/model1^2;
                        t22 = data2/model2^2;
                        for l = 0:NV-1
                            for m =l:NV-1
                                hessian(l*NV+m+1) = hessian(l*NV+m+1)+t21*newDudt1(l+1)*newDudt1(m+1)+t22*newDudt2(l+1)*newDudt2(m+1);
                                hessian(m*NV+l+1) = hessian(l*NV+m+1);
                            end
                        end
                        
                    end
                end
                %                                 dipshow(newDudt1);
            else
                mult = max( (1 + newLambda*10)/(1 + newLambda),1.3);
                newLambda = 10*newLambda;
%                 maxJump = maxJump*2;
                disp('CHOLERRER')
            end
            %                 kk
            
            
        end
        if kk == iterations
            error(tx) = newErr;
            num_of_iterations = kk;
        end
      %  plot_psf( PSF, newTheta(5), newTheta(6), newTheta(3), newTheta(4), phi0 )
    end
    
    P(tx, 1:length(newTheta)) = newTheta;
    P(tx, length(newTheta) + 1 : length(maxJump) + length(newTheta)) = maxJump;
    P(tx, length(maxJump) + length(newTheta) + 1) = num_of_iterations;
  
        
end









































