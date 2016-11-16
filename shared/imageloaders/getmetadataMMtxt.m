function mdo=getmetadataMMtxt(metafile,metadata)  
try
fid=fopen(metafile);
if fid>0
    minfo=fread(fid,[1,100000],'*char');
   calibfile='settings/CameraCalibration.xls';
    camcalib=readtable(calibfile);
    fclose(fid);
    md=minfoparsec(minfo,camcalib);
    mdo=md;
end

if nargin>2
    mdo=copyfields(metadata,mdo);
    fn=fieldnames(md);
    for k=1:length(fn)
        obj.metadata.assigned.(fn{k})=true;
    end
    mdo.allmetadata=copyfields(info,md);
    mdo.assigned.allmetadata=true;
%             obj.metadata.camerainfo=copyfields(obj.metadata.camerainfo,md);
            %determine if EM is used
            switch md.port
                case {'Conventional','Normal'}
                    mdo.EMon=false;
                    mdo.assigned.EMon=true;
                case {'Electron Multiplying', 'EM','Multiplication Gain'}
                    mdo.EMon=true;
                    mdo.EMon=true;
                otherwise 
                    md.port
                    mdo.EMon=true;      
            end
end
catch err
    imgs=dir([info.basefile filesep '*.tif']);
    mdo.allmetadata.files={imgs(:).name};
    mdo.allmetadata.frames=info.numberOfFrames;
    mdo.allmetadata.path=info.basefile;
    mdo.allmetadata.metafile='';
    disp('problems reading metadata')
end

end