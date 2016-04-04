function  [xfit,yfit]=LocUncFit(MLocUnc,MPConfInt,VminL,VminC);

%%%%%
%%%%%   This function fits the uncertainty curves. You may need to adjust
%%%%%   the following parameters for the best results: truncL, trincR,
%%%%%   hangout and the fiting parameters in cfobj. 
%%%%%

    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % %%%%%%% initialization parameters
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    truncL=0;      % number of data points to discard from begining of uncertainty data vs VminL
    truncR=5;      % number of data points to discard from end of uncertainty data vs VminL
    
    NC=length(VminC);  % number of minC values
    NL=length(VminL);  % number of minL values or equivalently length of VminL in the FOCALmain.m 
    hangout=10;    % data points in the range of  NL-truncR-hangout:NL-truncR  will be used for the estimation of asymptote and the SD in the asymptote.

    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % %%%%%%% fit the curves and plot them
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    for i=1:NC
           restline=mean(  MLocUnc(NL-truncR-hangout:NL-truncR,i)); % estimated asymptote  
           
           sminC=num2str(VminC(i)); 

           VLocUnc=(MLocUnc(:,i))';   
           VWeight=(MPConfInt(:,i))';
           
           truncVLocUnc=VLocUnc(1+truncL:end-truncR); 
           truncVWeight=VWeight(1+truncL:end-truncR);
           truncVminL=VminL(1+truncL:end-truncR);
           
           
                ftobj = fittype('Y0+A1.*exp(-1*((x-xL))./X0)');
                
                        % initial guess of fitting parameters
                        A1i=-((truncVLocUnc(8)-truncVLocUnc(1))/(truncVminL(8)-truncVminL(1))); 
                        xLi=2; 
                        X0i=4;                
                        y0i=min(truncVLocUnc);
                        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                
                cfobj = fit(truncVminL',truncVLocUnc',ftobj,'Weights', truncVWeight , 'Lower',[-100,1,4,1],'Upper',[100,20,20,3],'StartPoint',[A1i X0i y0i xLi])
           
          xfit=truncVminL(1):0.1:truncVminL(end);
          yfit=feval(cfobj,xfit); 
          
          % find the optimum minL
          Unc(i)=std(MLocUnc(end-truncR-hangout:end-truncR,i)); %  std of the last data points
          OminL(i)=log((  (Unc(i)+restline-cfobj.Y0) ./cfobj.A1)^-1)*cfobj.X0+cfobj.xL;
          %%%%%%%%%%%%%%%%%%%%%%%
          
          
          % calculate the uncertainty of the optimum minL        
          stdUnc(i) = Unc(i)*gamma((hangout-1)/2)/gamma(hangout/2)*(   (((hangout-1)/2)-(gamma(hangout/2)/gamma((hangout-1)/2))^2)    )^0.5; % std(std(of last data points))          
            UncUp(i)=(Unc(i)+restline-cfobj.Y0)+stdUnc(i);
            UncLow(i)=(Unc(i)+restline-cfobj.Y0)-stdUnc(i);          
          stdOminL(i)=max(  log((UncUp(i)./cfobj.A1)^-1)*cfobj.X0+cfobj.xL - OminL(i) ,  log((UncLow(i)./cfobj.A1)^-1)*cfobj.X0+cfobj.xL - OminL(i)   )   ;
          %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
      
          figure (1)
              if NC>1
                 subplot(3,ceil(NC/3),i)
              end
          hold on;  
              errorbar(VminL,VLocUnc,MPConfInt(:,i),'.k')
              plot(xfit,yfit,'r-'); 
              plot(xfit,(Unc(i)+restline),'g')   
                   title(['minC=' sminC])
                   xlabel('minL'); ylabel('Localization Uncertainty (nm)');    
                   xlim([VminL(1)-0.5 VminL(end)+0.5-truncR]);
                   sminLOpt = ['minL* = ' num2str(OminL(i))];
                   text(mean(VminL(1:end- truncR)),max(VLocUnc),sminLOpt)
          hold off;  
           
    end 

end


