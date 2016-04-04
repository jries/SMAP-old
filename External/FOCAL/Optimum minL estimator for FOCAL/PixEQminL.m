function [approvalList]=PixEQminL(FileName,PathName,minL)
%%%%%
%%%%%   This function output the index of those localizations that are
%%%%%   located in those pixels of density map that their value is equal to minL  
%%%%%

    global SRpixel Nneighbor XSRp YSRp;

    LocTable = importdata([PathName FileName]);

    Vx=LocTable.data(:,1); Vxind=max(0,floor( Vx/SRpixel))+1; % localization x coordinate in terms of SRpixels
    Vy=LocTable.data(:,2); Vyind=max(0,floor( Vy/SRpixel))+1; % localization y coordinate in terms of SRpixels

    [Nrows,Ncolumns]=size(LocTable.data); % row is the number of localizations; column is the number of saved parameters

    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % %%%%%%%  Construct the SRL Image 
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    SRLImage=zeros(YSRp,XSRp);

    for j=1:Nrows
        if Vyind(j)<=YSRp && Vxind(j)<=XSRp && Vyind(j)>0  && Vxind(j)>0
        SRLImage(Vyind(j),Vxind(j))=SRLImage(Vyind(j),Vxind(j))+1;   %   localization Image
        end
    end    
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % %%%%%%%  Construct the density map of SRL image:SRLDM
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    SRLDM= zeros(YSRp,XSRp);

    for i=2:YSRp-1
        for j=2:XSRp-1
            for m=-1:1
                 for n=-1:1
                     SRLDM(i,j)=SRLDM(i,j)+SRLImage(i+n,j+m);
                 end
            end
        end
    end    
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % %%%%%%%  find localizations within pixels with minL localizations
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ClusterMap= zeros(YSRp,XSRp);

    for i=1:Nrows;    
        if  SRLDM(Vyind(i),Vxind(i))==minL 
            approvalList(i)=1;
        else 
            approvalList(i)=0;
        end    
    end 
end    