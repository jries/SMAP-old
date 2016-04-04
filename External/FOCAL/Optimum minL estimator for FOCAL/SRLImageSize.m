function [ XSRp, YSRp ] = SRLImageSize( FileName,PathName )
%%%%%
%%%%%   This function extracts the size information of the original image 
%%%%%   from the header
%%%%%

    global pixel SRpixel;
    
    LocTable = importdata([PathName FileName]);

    StartIndex=strfind(LocTable.textdata{:},'max="');
    StopIndex=strfind(LocTable.textdata{:},'m" /');

    Xnm= str2num(    LocTable.textdata{1}(StartIndex(1)+5:StopIndex(1)-1)    ) *10^9; % X Size of camera image in terms of nm
    Ynm= str2num(    LocTable.textdata{1}(StartIndex(2)+5:StopIndex(2)-1)    ) *10^9; % Y Size of camera image in terms of nm

    Xp=round(1+ Xnm/pixel); % X Size of camera image in terms of pixels
    Yp=round(1+ Ynm/pixel); % Y Size of camera image in terms of pixels

    SRfactor=pixel/SRpixel;
    XSRp=floor(Xp*SRfactor);    % X Size of SR image in terms of SRpixels
    YSRp=floor(Yp*SRfactor);    % Y Size of SR image in terms of SRpixels

end

