function [BPathName, BFileName] = ConfineLocTable( FileName,PathName )
%%%%%
%%%%%    This function imports the localization table and filters out
%%%%%    localizations that are out of the image boundaries. 
%%%%%    A localization can fall outside of the boundaries of an image in 
%%%%%    the process of drift correction 
%%%%%

    global SRpixel XSRp YSRp;

    LocTable = importdata([PathName FileName]);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%    Find pixel index of localizations
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    Vx=LocTable.data(:,1); Vxind=max(0,floor( Vx/SRpixel))+1;   % localization x coordinate in terms of SRpixels
    Vy=LocTable.data(:,2); Vyind=max(0,floor( Vy/SRpixel))+1;   % localization y coordinate in terms of SRpixels

    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % %%%%%%%  Reject the out of bound localizations and save the rest in a new
    % %%%%%%%  file
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    mkdir(PathName,FileName(1:end-4));
    BFileName='BLocTable.txt';
    BPathName=fullfile(PathName,FileName(1:end-4),filesep);
    Bfile=fullfile(BPathName,BFileName);
    fileid=fopen(Bfile,'w');
    fprintf(fileid,'%s\r\n',LocTable.textdata{1}(:));   
    [Nrows,Ncolumns]=size(LocTable.data);

    for i=1:Nrows
        if Vyind(i)<=YSRp && Vxind(i)<=XSRp && Vyind(i)>0  && Vxind(i)>0
        Linei=[];
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

