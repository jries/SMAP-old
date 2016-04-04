function  [FOCALPathName,FOCALFileName]=LocToROI(FileName,PathName,minC);

%%%%%
%%%%%   This function uses the density filtered localization table,
%%%%%   identifies clusters, rejects small clusters (< minC (10 in FOCAL))
%%%%%   and saves a new localization file from the localizations within
%%%%%   accepted clusters.   
%%%%%   Number of clusters and the index of localizations (NROIs,ROIsIndList)can be used in the output argument whenever needed. In this caseThe FOCALmain.m should be changed accordingly.

    global SRpixel Nneighbor XSRp YSRp;
    localizations = [PathName FileName];
    LocTable = importdata(localizations);

    Vx=LocTable.data(:,1); Vxind=max(0,floor( Vx/SRpixel))+1; % localization x coordinate in terms of SRpixels
    Vy=LocTable.data(:,2); Vyind=max(0,floor( Vy/SRpixel))+1; % localization y coordinate in terms of SRpixels
    [Nrows,Ncolumns]=size(LocTable.data); % row is the number of localizations; column is the number of saved parameters

    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % %%%%%%%  Construct the Black and White super resolution Image
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        BWSRImage=zeros(YSRp,XSRp);

        for j=1:Nrows        
            BWSRImage(Vyind(j),Vxind(j))=255;        
        end

    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % %%%%%%%  find clusters(ROIs) and make their list
    % %%%%%%%  
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    bwSRImage = im2bw(BWSRImage,0.1);
    if Nneighbor~=4
       Nneighbor=8;
    end

    ccSRImage = bwconncomp(bwSRImage,Nneighbor);
    ROIsIndList=[];
    N=0;
    for i=1:length(ccSRImage.PixelIdxList)            
                               
        if length(ccSRImage.PixelIdxList{i})>=minC % rejects clusters that are smaller than minC pixels
            [subx,suby]= ind2sub([YSRp,XSRp],ccSRImage.PixelIdxList{i});                        
            if min(subx)~=max(subx) && min(suby)~=max(suby)  % Excludes clusters that are 1D and not 2D
                N=N+1; 
                ROIsIndList{N}=ccSRImage.PixelIdxList{i};         
            end
        end
    end

    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % %%%%%%% Cunstruct a mask Image from the accepted ROIs. The mask will
    % %%%%%%% be used to identify localizations that are
    % %%%%%%% member of any of the clusters and to Generate an approval
    % %%%%%%% list. Other localizations are noise and will be rejected.
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    ROIImage=zeros(YSRp,XSRp);
    for i=1:length(ROIsIndList)
        ROIImage(ROIsIndList{i})=255;
    end                                

    for i=1:Nrows;
        approval(i)=0;
        if ROIImage(Vyind(i),Vxind(i))==255   
           approval(i)=1;
        end
    end

    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % %%%%%%% Save the clustered localization table by FOCAL
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    sminC=num2str(minC);

    FOCALFileName=['FOCAL_minC' sminC '_' FileName];
    FOCALPathName=PathName;
    FOCALfile=fullfile(FOCALPathName,FOCALFileName);

    fileid=fopen(FOCALfile,'w');
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