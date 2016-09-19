function mdo=getmetadataMM(obj)
info=getimageinfo(obj.file);
            
            fid=fopen(info.metafile);
            if fid>0
                minfo=fread(fid,[1,100000],'*char');
                camcalib=readtable(obj.calfile);
                fclose(fid);
                md=minfoparsec(minfo,camcalib);
            end
             obj.metadata.numberOfFrames=info.numberOfFrames;
             obj.metadata.basefile=info.basefile;
            obj.metadata=copyfields(obj.metadata,md);
            fn=fieldnames(md);
            for k=1:length(fn)
                obj.metadata.assigned.(fn{k})=true;
            end
            obj.metadata.allmetadata=copyfields(info,md);
             obj.metadata.assigned.allmetadata=true;
%             obj.metadata.camerainfo=copyfields(obj.metadata.camerainfo,md);
            %determine if EM is used
            switch md.port
                case {'Conventional','Normal'}
                    obj.metadata.EMon=false;
                    obj.metadata.assigned.EMon=true;
                case {'Electron Multiplying', 'EM','Multiplication Gain'}
                    obj.metadata.EMon=true;
                    obj.metadata.assigned.EMon=true;
                otherwise 
                    md.port
                    obj.metadata.EMon=true;      
            end
            mdo=obj.metadata;
end