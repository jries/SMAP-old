function filestruct=locsaveFileinfo(obj)
tifs=struct('image',[],'info',[]);
finfo=obj.getPar('loc_fileinfo');
sfile=finfo.basefile;
filename=[sfile '_sml.mat'];  
infost=obj.getPar('loc_cameraSettings');
infost=copyfields(infost,finfo);

%git version
[gitstatus,head]=system('git rev-parse HEAD');
if gitstatus == 0
infost.git.head=head;
else
    infost.git.head=[];
end
[gitstatus,branch]=system('git branch');
if gitstatus== 0
ind=strfind(branch,'*');
le=find(branch(ind+1:end)==10,1,'first')+ind;
infost.git.branch=branch(ind+2:le-1);
else
infost.git.branch='not found';
end

filestruct=struct('info',infost,'average',[],'name',filename,...
    'number',obj.filenumber,...
    'numberOfTif',0,'tif',tifs,'raw',struct('image',[],'frame',[]));
if ~isfield(filestruct.info,'roi')||isempty(filestruct.info.roi)
    try
        mm=imfinfo(finfo.filename);
        obj.fileinfo.roi=[0 0 mm(1).Width mm(1).Height];

    catch
        try
            filestruct.info.roi=[0 0 filestruct.info.Width filestruct.info.Height];
        catch
            filestruct.info.roi=[0 0 512 512];
        end
    end
    filestruct.info.roi=obj.fileinfo.roi;
end
