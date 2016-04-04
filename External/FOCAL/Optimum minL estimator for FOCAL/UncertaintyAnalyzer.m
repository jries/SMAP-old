function [LocUnc,PCI]=UncertaintyAnalyzer(FileName,PathName,approvalList);

%%%%%
%%%%%   This function calculates average localization uncertainty by the 
%%%%%   method described in "Endesfelder, Ulrike, et al. "A simple method 
%%%%%   to estimate the average localization precision of a single-molecule
%%%%%   localization microscopy experiment." Histochemistry and cell 
%%%%%   biology 141.6 (2014): 629-638.".   
%%%%%

    global pixel;
    LocTable = importdata([PathName FileName]);

    Vx=LocTable.data(:,1); 
    Vy=LocTable.data(:,2); 
    frames=LocTable.data(:,3);

    [Nrows,Ncolumns]=size(LocTable.data); % row is the number of localizations; column is the number of saved parameters

    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % %%%%%%% Find pairs of consequetive overlaped blinks and their mutual
    % %%%%%%% distance  
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    Dframes=diff(frames);
    transitionInd=find(Dframes>0);  % find where in the localization table the frame number changes
    transitionInd=[0 transitionInd' Nrows];

    VSTAdx=[];    VSTAdy=[]; % Arrays of dx and dy of Spatial and Temporal Adjacent localizations

    for i=2:length(transitionInd)-1
        if Dframes(transitionInd(i)) ==1 % if the frame number increments by 1 (temporal neighbors) find adjacent spatial neighbor localizations
            Ap=approvalList(transitionInd(i-1)+1 : transitionInd(i));             
            pVx=Vx(transitionInd(i-1)+1 : transitionInd(i)); % p stands for previous. pVx is 
            pVy=Vy(transitionInd(i-1)+1 : transitionInd(i));
            for j=transitionInd(i)+1:transitionInd(i+1)
                VDxy=(  ((pVx-Vx(j)).^2+(pVy-Vy(j)).^2).^0.5  ).* (max(Ap,approvalList(j)))'; %if one of the former or later localizations is in the approval list: VDxy is calculated otherwise it will be set to 0
                SpatialAdjacentLocs=find(VDxy<pixel/2 & VDxy>0);
                if length(SpatialAdjacentLocs)==1
                    VSTAdx=[VSTAdx pVx(SpatialAdjacentLocs)-Vx(j)];
                    VSTAdy=[VSTAdy pVy(SpatialAdjacentLocs)-Vy(j)];
                end
            end
        end

    end
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % %%%%%%% Calculate the average localization uncertainty
    % %%%%%%%                  
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
    N=length(VSTAdx);
    if N>=1
        Vd2=VSTAdx.^2+VSTAdy.^2;
        Vd1=Vd2.^0.5;
        custpdf = @(dij,sig) (dij/(2*sig^2)).*exp(-1*dij.^2/(4*sig^2));
        [LocUnc, PCI2] = mle(Vd1,'pdf',custpdf,'start',15);
        PCI=(PCI2(2)-PCI2(1))/2;   
    else
        LocUnc=0;
        PCI=0;
    end    
    
end
