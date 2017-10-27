function [P,CRLB, LL,update, error] =  kernel_MLEfit_Spline_LM_multichannel_finalized(d_data,coeff, shared,dTAll,iterations)
% update computeDelta3Dj and fAt3Dj. Add kernel_derivativeSpline
% make it capable of handling coeff with random size using Spline3D_v2
%work with odd number of ROI for calibration
%try to solve the negative value problem


% maxJump_Init = [1, 1,spline_zsize/5,100,100];


pi = single(3.141592);
spline_xsize = size(coeff, 1);
spline_ysize = size(coeff, 2);
spline_zsize = size(coeff, 3);
PSFSigma = single(1.5);
noChannels = size(d_data, 4);
sz= size(d_data,1);

NV = 5*noChannels-sum(shared(:))*(noChannels-1);

maxJump_Init = zeros(5,noChannels);
for i = 1:noChannels
    maxJump_Init(:,i)=[1; 1;spline_zsize/5;100;100];
    coeffall{i} = coeff(:,:,:,:,i);
end


delta_f1 = zeros(64,1);
delta_dxf1 = zeros(64,1);
delta_dyf1 = zeros(64,1);
delta_dzf1 = zeros(64,1);


% n=1;
% maxJump = zeros(NV,1);
% for i = 1:5
%     if shared(i)==1
%         maxJump(n)=maxJump_Init(i);
%     else
%         for j= 1:noChannels
%             maxJump(n+j-1)=maxJump_Init(i);
%         end
%        n = n+j-1;
%     end
%     n=n+1;
% end

Nfits = size(d_data,3);

xc = zeros(1,noChannels);
yc = zeros(1,noChannels);
zc = zeros(1,noChannels);


delta_f = zeros(64,noChannels);
delta_dxf = zeros(64,noChannels);
delta_dyf = zeros(64,noChannels);
delta_dzf = zeros(64,noChannels);

model = zeros(1,noChannels);
data = zeros(1,noChannels);
newDudt = zeros(5,noChannels);


% newDudtall = zeros(NV,noChannels);


% dudt = single(zeros(sz,sz,NV));

newTheta = single(zeros(5,noChannels));

newThetaAll= zeros(NV,1);
% newTheta2 = single(zeros(5,1));

% oldTheta = newTheta;

% M = zeros(NV*NV,1);
% Minv = zeros(NV*NV,1);

Nfits = size(d_data,3);

error = zeros(Nfits, 1);
model1 =0;
newDudt1=0;
tolerance=1e-6;

num_of_iterations = 0;

% sigma_x = zeros(Nfits, 1);
% sigma_y = zeros(Nfits, 1);
% sigma_to_z = 0.55;

% weightPara = [1,1;1,1;1,1;1,0;0,1;1,0;0,1];
for tx = 1:Nfits
    tx;
    newLambda = 0.01;
    %newpeak
    dT = dTAll(:,:,tx);
    for i = 1:noChannels
        [newTheta(1,i),newTheta(2,i)]=kernel_CenterofMass2D(sz,single(d_data(sz*sz*(tx-1)+(i-1)*sz*sz*Nfits+1:sz*sz*(tx)+(i-1)*sz*sz*Nfits)));
        
        
        [Nmax,newTheta(5,i)] = kernel_GaussFMaxMin2D(sz,PSFSigma,single(d_data(sz*sz*(tx-1)+(i-1)*sz*sz*Nfits+1:sz*sz*(tx)+(i-1)*sz*sz*Nfits)));
        
        newTheta(4,i)=max(0,(Nmax-newTheta(5,i))*2*pi*PSFSigma*PSFSigma);
        %         newTheta(5)=max(0,(Nmax2-newTheta(7))*2*pi*PSFSigma*PSFSigma);
        
        newTheta(3,i) = spline_zsize/2;
        
        maxJump_Init(4,i) =  max(maxJump_Init(4,i),newTheta(4,i));
        maxJump_Init(5,i) =  max(maxJump_Init(5,i),newTheta(5,i));
    end
    %     newTheta2(4) =  newTheta1(4);
    %     newTheta(5) = in_pos;
    %     newTheta(6) = 0;
    
    %     newLambda = (1);
    

    
    for i = 1:5
        if shared(i)==1
            for j = 2:noChannels
                newTheta(i,j)= newTheta(i,1)+dT(i,j);
            end
        end
    end
    
    newUpdate = 100000*ones(NV,1);
    oldUpdate = 100000*ones(NV,1);
    
    %     maxJump = [1, 1,spline_zsize/5,100,100, 100,100];
    %     maxJump(4) = max(maxJump(4),newTheta(4));
    %     maxJump(5) = max(maxJump(5),newTheta(5));
    %
    %     maxJump(6) = max(maxJump(6),newTheta(6));
    %     maxJump(7) = max(maxJump(7),newTheta(7));
    
    n=1;
    maxJump = zeros(NV,1);
    for i = 1:5
        if shared(i)==1
            maxJump(n)=mean(maxJump_Init(i,:));
            newThetaAll(n)=newTheta(i,1);
        else
            for j= 1:noChannels
                maxJump(n+j-1)=maxJump_Init(i,j);
                newThetaAll(n+j-1)=newTheta(i,j);
            end
            n = n+j-1;
        end
        n=n+1;
    end

    
    
    %updateFitValues3D
%     xc = single(-2*(newTheta(1) - sz/2+0.5));
%     yc = single(-2*(newTheta(2) - sz/2+0.5));


    off = floor(((spline_xsize+1)-sz)/2);
    for i = 1:noChannels
        xc(i)=single(-1 * (newTheta(1,i) - sz / 2 + 0.5));
        yc(i) = single(-1 * (newTheta(2,i) - sz / 2 + 0.5));
        zc(i) = single(newTheta(3,i)-floor(newTheta(3,i)));
        
        xstart(i)=floor(xc(i));
        ystart(i)=floor(yc(i));
        zstart(i)=floor(newTheta(3,i));
        
        xc(i) = xc(i)-floor(xc(i));
        yc(i) = yc(i)-floor(yc(i));
        
    end
    
    jacobian = single(zeros(NV, 1));
    hessian = single(zeros(NV, NV));
    
    for i = 1: noChannels
        [delta_f(:,i),delta_dxf(:,i),delta_ddxf1,delta_dyf(:,i),delta_ddyf1,delta_dzf(:,i),delta_ddzf1] = computeDelta3Dj_v2(single(xc(i)),single(yc(i)),single(zc(i)));
        
        %      [delta_f2,delta_dxf2,delta_ddxf2,delta_dyf2,delta_ddyf2,delta_dzf2,delta_ddzf2] = computeDelta3Dj_v2(single(xc2),single(yc2),single(zc));
    end
    
%      newTheta1 = [newTheta(1),newTheta(2),newTheta(3),newTheta(4),newTheta(6)];
%      newTheta2 = [newThetaT(1),newThetaT(2),newTheta(3),newTheta(5),newTheta(7)];
    newDudtAll = zeros(NV,noChannels); 
    newErr = 0;
    for ii = 0:sz-1
        for jj = 0:sz-1
            
            
            for i = 1:noChannels
                data(i) = single(d_data(sz*sz*(tx-1)+sz*jj+ii+1+sz*sz*Nfits*(i-1)));
                %                 data2 = single(d_data2(sz*sz*(tx-1)+sz*jj+ii+1));
           
                
                delta_f1 = delta_f(:,i);
                delta_dxf1 = delta_dxf(:,i);
                delta_dyf1 = delta_dyf(:,i);
                delta_dzf1 = delta_dzf(:,i);
%                 coeff1 = coeff(:,:,:,:,i);
                newTheta1 = newTheta(:,i);
                
                [newDudt1, model1] =  kernel_DerivativeSpline_v2_finalized(ii+xstart(i)+off,jj+ystart(i)+off,zstart(i),spline_xsize,spline_ysize,spline_zsize,delta_f1,delta_dxf1,delta_dyf1,delta_dzf1,coeffall{i},newTheta1,NV);
                newDudt(:,i)=newDudt1;
                model(i)=model1;
                %                 [newDudt2, model2] =  kernel_DerivativeSpline_v2(ii+xstart2+off,jj+ystart2+off,zstart,spline_xsize,spline_ysize,spline_zsize,delta_f2,delta_dxf2,delta_dyf2,delta_dzf2,coeff2,newTheta2,2);
            end
            
            
            n=1;
            for i = 1:5
                if shared(i)==1
                    newDudtAll(n,:)=newDudt(i,:);
                else
                    for j= 1:noChannels
                        newDudtAll(n+j-1,j)=newDudt(i,j);
                    end
                    n = n+j-1;
                end
                n=n+1;
            end
            
            %             if data>0
            %                 newErr = newErr +2*((model-data)-data*log(model/data));
            %             elseif data ==0
            %                 newErr = newErr + 2*model;
            %             end
            for i = 1:noChannels
                if data(i)>0
                    newErr = newErr +2*((model(i)-data(i))-data(i)*log(model(i)/data(i)));
                else
                    newErr = newErr + 2*model(i);
                    data = 0;
                end
            end
            
            
            
            
            t1 = single(1 - data./ model);
            
            for l = 1:NV
                for j= 1:noChannels
                    jacobian(l) = jacobian(l)+t1(j)*newDudtAll(l,j);
                end
                %                 jacobian2(l) = jacobian2(l)+t12*newDudt2(l);
            end
            
             
            t2 = data./model.^2;
            for l = 0:NV-1
                for m =l:NV-1
                    for j = 1:noChannels
                        hessian(l*NV+m+1) = hessian(l*NV+m+1)+t2(j)*newDudtAll(l+1,j)*newDudtAll(m+1,j);
                    end
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
        
        tx
        kk
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
                    newThetaAll(ll) = newThetaAll(ll)-newUpdate(ll);
                    update(kk,ll) = newUpdate(ll);
                end
                
                
                n=1;
                for i = 1:5
                    if shared(i)==1
                        %                         maxJump(n)=mean(maxJump_Init(i,:));
                        switch i
                            case 1
                            case 2
                                newThetaAll(n)= max(newThetaAll(n),(sz-1)/2-sz/4.0);
                                newThetaAll(n) = min(newThetaAll(n),(sz-1)/2+sz/4.0);
                            case 3
                                newThetaAll(n) = max(newThetaAll(n), 0);
                                newThetaAll(n) = min(newThetaAll(n), spline_zsize);
                            case 4
                                newThetaAll(n) = max(newThetaAll(n), 1);
                            case 5
                                newThetaAll(n) = max(newThetaAll(n), 0.01);
                        end
                        
                        for j = 1:noChannels
                            newTheta(i,j)=newThetaAll(n)+dT(i,j);
                        end
                    else
                        for j= 1:noChannels
                            %                             newThetaAll(n+j-1)=newTheta(i,1);
                            switch i
                                case 1
                                case 2
                                    newThetaAll(n+j-1)= max(newThetaAll(n+j-1),(sz-1)/2-sz/4.0);
                                    newThetaAll(n+j-1) = min(newThetaAll(n+j-1),(sz-1)/2+sz/4.0);
                                case 3
                                    newThetaAll(n+j-1) = max(newThetaAll(n+j-1), 0);
                                    newThetaAll(n+j-1) = min(newThetaAll(n+j-1), spline_zsize);
                                case 4
                                    newThetaAll(n+j-1) = max(newThetaAll(n+j-1), 1);
                                case 5
                                    newThetaAll(n+j-1) = max(newThetaAll(n+j-1), 0.01);
                            end
                            newTheta(i,j)=newThetaAll(n+j-1);
                            
                        end
                        n = n+j-1;
                    end
                    n=n+1;
                end
                
                
                
                %                 newTheta(1) = max(newTheta(1),(sz-1)/2-sz/4.0);
                %                 newTheta(1) = min(newTheta(1),(sz-1)/2+sz/4.0);
                %                 newTheta(2) = max(newTheta(2),(sz-1)/2-sz/4.0);
                %                 newTheta(2) = min(newTheta(2),(sz-1)/2+sz/4.0);
                %                 newTheta(4) = max(newTheta(4), 1);
                %                 newTheta(5) = max(newTheta(5), 1);
                %                 newTheta(6) = max(newTheta(6), 0.01);
                %                 newTheta(7) = max(newTheta(7), 0.01);
                %                 newTheta(3) = max(newTheta(3), 0);
                %                 newTheta(3) = min(newTheta(3), spline_zsize);
                
                
                
                for i = 1:noChannels
                    xc(i)=single(-1 * (newTheta(1,i) - sz / 2 + 0.5));
                    yc(i) = single(-1 * (newTheta(2,i) - sz / 2 + 0.5));
                    zc(i) = single(newTheta(3,i)-floor(newTheta(3,i)));
                    
                    xstart(i)=floor(xc(i));
                    ystart(i)=floor(yc(i));
                    zstart(i)=floor(newTheta(3,i));
                    
                    xc(i) = xc(i)-floor(xc(i));
                    yc(i) = yc(i)-floor(yc(i));
                    
                end
                
                jacobian = single(zeros(NV, 1));
                hessian = single(zeros(NV, NV));
                
                for i = 1: noChannels
                    [delta_f(:,i),delta_dxf(:,i),delta_ddxf1,delta_dyf(:,i),delta_ddyf1,delta_dzf(:,i),delta_ddzf1] = computeDelta3Dj_v2(single(xc(i)),single(yc(i)),single(zc(i)));
                    
                    %      [delta_f2,delta_dxf2,delta_ddxf2,delta_dyf2,delta_ddyf2,delta_dzf2,delta_ddzf2] = computeDelta3Dj_v2(single(xc2),single(yc2),single(zc));
                end
                
                %      newTheta1 = [newTheta(1),newTheta(2),newTheta(3),newTheta(4),newTheta(6)];
                %      newTheta2 = [newThetaT(1),newThetaT(2),newTheta(3),newTheta(5),newTheta(7)];
                newDudtAll = zeros(NV,noChannels);
                newErr =0;
                for ii = 0:sz-1
                    for jj = 0:sz-1
                        
                        
                        for i = 1:noChannels
                            data(i) = single(d_data(sz*sz*(tx-1)+sz*jj+ii+1+sz*sz*Nfits*(i-1)));
                            %                 data2 = single(d_data2(sz*sz*(tx-1)+sz*jj+ii+1));
                            
                            delta_f1 = delta_f(:,i);
                            delta_dxf1 = delta_dxf(:,i);
                            delta_dyf1 = delta_dyf(:,i);
                            delta_dzf1 = delta_dzf(:,i);
%                             coeff1 = coeff(:,:,:,:,i);
                            newTheta1 = newTheta(:,i);
                            
                            [newDudt1, model1] =  kernel_DerivativeSpline_v2_finalized(ii+xstart(i)+off,jj+ystart(i)+off,zstart(i),spline_xsize,spline_ysize,spline_zsize,delta_f1,delta_dxf1,delta_dyf1,delta_dzf1,coeffall{i},newTheta1,NV);
                            newDudt(:,i)=newDudt1;
                            model(i)=model1;
                            %[newDudt(:,i), model(i)] =  kernel_DerivativeSpline_v2_finalized(ii+xstart(i)+off,jj+ystart(i)+off,zstart(i),spline_xsize,spline_ysize,spline_zsize,delta_f(:,i),delta_dxf(:,i),delta_dyf(:,i),delta_dzf(:,i),coeff(:,:,:,:,i),newTheta(:,i),NV);
                            
                          end
                        
                        
                        n=1;
                        for i = 1:5
                            if shared(i)==1
                                newDudtAll(n,:)=newDudt(i,:);
                            else
                                for j= 1:noChannels
                                    newDudtAll(n+j-1,j)=newDudt(i,j);
                                end
                                n = n+j-1;
                            end
                            n=n+1;
                        end
                        
                        %             if data>0
                        %                 newErr = newErr +2*((model-data)-data*log(model/data));
                        %             elseif data ==0
                        %                 newErr = newErr + 2*model;
                        %             end
                        for i = 1:noChannels
                            if data(i)>0
                                newErr = newErr +2*((model(i)-data(i))-data(i)*log(model(i)/data(i)));
                            else
                                newErr = newErr + 2*model(i);
                                data = 0;
                            end
                        end
                        
                        
                        
                        
                        t1 = single(1 - data./ model);
                        
                        for l = 1:NV
                            for j= 1:noChannels
                                jacobian(l) = jacobian(l)+t1(j)*newDudtAll(l,j);
                            end
                            %                 jacobian2(l) = jacobian2(l)+t12*newDudt2(l);
                        end
                        
                        
                        t2 = data./model.^2;
                        for l = 0:NV-1
                            for m =l:NV-1
                                for j = 1:noChannels
                                    hessian(l*NV+m+1) = hessian(l*NV+m+1)+t2(j)*newDudtAll(l+1,j)*newDudtAll(m+1,j);
                                end
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
    
    
    for i = 1:noChannels
        xc(i)=single(-1 * (newTheta(1,i) - sz / 2 + 0.5));
        yc(i) = single(-1 * (newTheta(2,i) - sz / 2 + 0.5));
        zc(i) = single(newTheta(3,i)-floor(newTheta(3,i)));
        
        xstart(i)=floor(xc(i));
        ystart(i)=floor(yc(i));
        zstart(i)=floor(newTheta(3,i));
        
        xc(i) = xc(i)-floor(xc(i));
        yc(i) = yc(i)-floor(yc(i));
        
    end
    
    for i = 1: noChannels
        [delta_f(:,i),delta_dxf(:,i),delta_ddxf1,delta_dyf(:,i),delta_ddyf1,delta_dzf(:,i),delta_ddzf1] = computeDelta3Dj_v2(single(xc(i)),single(yc(i)),single(zc(i)));
        
        %      [delta_f2,delta_dxf2,delta_ddxf2,delta_dyf2,delta_ddyf2,delta_dzf2,delta_ddzf2] = computeDelta3Dj_v2(single(xc2),single(yc2),single(zc));
    end
    
    M = zeros(NV,NV);
    Div = 0;
    for ii = 0:sz-1
        for jj = 0:sz-1
            
            
            for i = 1:noChannels
                data(i) = single(d_data(sz*sz*(tx-1)+sz*jj+ii+1+sz*sz*Nfits*(i-1)));
                %                 data2 = single(d_data2(sz*sz*(tx-1)+sz*jj+ii+1));
                
                delta_f1 = delta_f(:,i);
                delta_dxf1 = delta_dxf(:,i);
                delta_dyf1 = delta_dyf(:,i);
                delta_dzf1 = delta_dzf(:,i);
%                 coeff1 = coeff(:,:,:,:,i);
                newTheta1 = newTheta(:,i);
                
                [newDudt1, model1] =  kernel_DerivativeSpline_v2_finalized(ii+xstart(i)+off,jj+ystart(i)+off,zstart(i),spline_xsize,spline_ysize,spline_zsize,delta_f1,delta_dxf1,delta_dyf1,delta_dzf1,coeffall{i},newTheta1,NV);
                newDudt(:,i)=newDudt1;
                model(i)=model1;
                %[newDudt(:,i), model(i)] =  kernel_DerivativeSpline_v2_finalized(ii+xstart(i)+off,jj+ystart(i)+off,zstart(i),spline_xsize,spline_ysize,spline_zsize,delta_f(:,i),delta_dxf(:,i),delta_dyf(:,i),delta_dzf(:,i),coeff(:,:,:,:,i),newTheta(:,i),NV);
                
                %                 [newDudt2, model2] =  kernel_DerivativeSpline_v2(ii+xstart2+off,jj+ystart2+off,zstart,spline_xsize,spline_ysize,spline_zsize,delta_f2,delta_dxf2,delta_dyf2,delta_dzf2,coeff2,newTheta2,2);
            end
            
            
            n=1;
            for i = 1:5
                if shared(i)==1
                    newDudtAll(n,:)=newDudt(i,:);
                else
                    for j= 1:noChannels
                        newDudtAll(n+j-1,j)=newDudt(i,j);
                    end
                    n = n+j-1;
                end
                n=n+1;
            end
            
            
            for l = 0:NV-1
                for m =l:NV-1
                    for j = 1:noChannels
                        M(l*NV+m+1) = M(l*NV+m+1)+newDudtAll(l+1,j)*newDudtAll(m+1,j)/model(j);
                    end
                    M(m*NV+l+1) = M(l*NV+m+1);
                end
            end
            
            
            for i = 1:noChannels
                if data(i)>0
                    Div = Div +data(i)*log(model(i))-model(i)-data(i)*log(data(i))+data(i);
                else
                    Div = Div -model(i);
%                     data = 0;
                end
            end                    
        end
    end
    
    Minv = inv(M);
    if ~isempty(kk)
        iteration = kk;
    else
        iteration = 0;
    end
    for kk = 1:NV
        P(Nfits*(kk-1)+tx) = newThetaAll(kk);
        P(Nfits*(NV+kk)+tx) = maxJump(kk);
        CRLB(Nfits*(kk-1)+tx)=Minv(kk,kk);
    end
    LL(tx,1) = double(Div);
    LL(tx,2) = chi2pdf(-2*Div,sz*sz-(5));
    P(Nfits*(NV)+tx) = iteration;
   
  
        
end

P = reshape(P,Nfits,NV+1+NV);
CRLB = reshape(CRLB,Nfits,NV);







































