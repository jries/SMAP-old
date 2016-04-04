function [DFPathName,DFFileName]=DensityFilter(FileName,PathName,minL)
%%%%%
%%%%%   This function implements a density filter on the Localization Table
%%%%%   and saves the new table. Localizations within Core and Border
%%%%%   pixels will be approved and the remining localization will be 
%%%%%   rejected as noise 
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
    % %%%%%%%  construct a map of Core pixels in Clusters: Cores=1, others=0
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    CoreMap= zeros(YSRp,XSRp);

    for i=1:Nrows;    
        if  SRLDM(Vyind(i),Vxind(i))>=minL
            CoreMap(Vyind(i),Vxind(i))=1;
        end    
    end
       
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % %%%%%%%  1. Find localizations in the border and Core pixels and add them to
    % %%%%%%%  the approval list.
    % %%%%%%%  2. Enhance the Cluster map by adding borders. Cores and borders=1, others=0
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    if Nneighbor~=4 
       Nneighbor=8; 
    end

    CN=[0,0;-1,0;0,-1;1,0;0,1;-1,-1;1,1;-1,1;1,-1]; % relative coordinate of Connected Neighbors

    for i=1:Nrows;
        approval(i)=0;
            for j=1:Nneighbor+1;
                   if CoreMap(Vyind(i)+CN(j,1),Vxind(i)+CN(j,2))>0   
                      approval(i)=1;                        
                   end
            end
    end

    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % %%%%%%%  save the density filtered localization table
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    sminL=num2str(minL); 

    DFFileName=['DF' '_minL' sminL '.txt' ];
    DFPathName=PathName;
    DFfile=fullfile(DFPathName,DFFileName);

    fileid=fopen(DFfile,'w');
    fprintf(fileid,'%s\r\n',LocTable.textdata{1}(:));   

    for i=1:Nrows;
        Linei=[];
            if approval(i)==1
                 for k=1:Ncolumns-1 
                     C{k}=num2str(LocTable.data(i,k));
                     Linei=[Linei,C{k},' '];
                 end
                    k=Ncolumns;
                    C{k}=num2str(LocTable.data(i,k));
                    Linei=[Linei,C{k}];
                    fprintf(fileid,'%s\r\n',Linei);                
            end        
    end
    fclose('all');

end    