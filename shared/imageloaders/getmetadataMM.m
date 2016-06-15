function mdo=getmetadataMM(obj)
info=getimageinfo(obj.file);
            
            fid=fopen(info.metafile);
            if fid>0
                minfo=fread(fid,[1,100000],'*char');
                camcalib=readtable(obj.calfile);
                fclose(fid);
                md=minfoparsec(minfo,camcalib);
            end
            obj.metadata.frames=info.numberOfFrames;
            obj.metadata=copyfields(obj.metadata,md);
            obj.metadata.allmetadata=copyfields(info,md);
            obj.metadata.camerainfo=copyfields(obj.metadata.camerainfo,md);
            %determine if EM is used
            switch md.port
                case {'Conventional','Normal'}
                    obj.metadata.em=false;
                case {'Electron Multiplying', 'EM'}
                    obj.metadata.em=true;
                otherwise 
                    md.port
                    obj.metadata.em=true;      
            end
            mdo=obj.metadata;
end