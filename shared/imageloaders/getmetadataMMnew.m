function mdo=getmetadataMMnew(metafile)
            fid=fopen(metafile);
            if fid<=0
                mdo=[];
                return
            end
            
            
            minfo=fread(fid,[1,100000],'*char');
            nread=100000;
            fseek(fid,-nread,'eof');
            minfo2=fread(fid,[1,nread],'*char');
            fclose(fid);
            tt=textscan(minfo,'%s','delimiter',',');
            t=tt{1};
            md=cell(length(t),2);
            indg=false(length(t),1);
            for k=1:length(t)
                x=textscan(t{k},'%s','delimiter',':' );
                xh=x{1};
                if length(xh)>1
                    indg(k)=true;
                    md(k,1)=strrep(xh(1),'"','');
                    md(k,2)=strrep(xh(2),'"','');
                end
                
            end
            md=md(indg,:);
            [~,ic]=unique(md(:,1));
            mdo=md(ic,:);   
            %remove frame keys
            mdo(strncmp(mdo(:,1),'FrameKey',8),:)=[];
            
            %manual
             ind=strfindfast(minfo,'"ROI": [',1);
             troi=textscan(minfo(ind+10:ind+100),'%d','delimiter',',');
             mdo(end+1,:)={'ROI direct',num2str(troi{:}')};

             
            ind=strfindfast(minfo2,'"FrameKey-',1,-1);

            ind2=strfindfast(minfo2,'"',ind);
            str=minfo2(ind:ind2-1);
            ind3=strfind(str,'-');
            str2=str(ind3(end)+1:end);
            numberOfFrames=str2double(str2)+1;
            mdo(end+1,:)={'numberOfFrames direct',num2str(numberOfFrames)};
             %sort
             [~,inds]=sortrows(mdo,1);
             mdo=mdo(inds,:);
end