classdef metadataSMAP<handle
    %metadata class for smap. mirrors the locData.files.file.info
    
    properties
        Width=512;
        Height=512;
        roi=[];
        camerainfo=struct('camID','','preAmp','','port','','temperature',0,'readoutrate','');
        allmetadata=[];
        exposure=1;
        emgain=1;
        EMon=false;
        conversion=1;
        offset=100;
        pixsize=0.1;
        timediff=30;
        comment='';
        numberOfFrames=0;
        pix2phot=[];
        basefile;
    end
    
    methods
        function roi=get.roi(obj)
            if isempty(obj.roi)
                roi=[0 0 obj.Width obj.Height];
            else
                roi=obj.roi;
            end
        end
        function cv=get.pix2phot(obj)
            if isempty(obj.pix2phot)
                if obj.EMon
                    cv=obj.conversion/obj.emgain;
                else
                    cv=obj.conversion;
                end
                obj.pix2phot=cv;
            else
                cv=obj.pix2phot;
            end
        end
        function set.roi(obj,roi)
            if isempty(roi)
                return
            end
            obj.roi=roi;
            obj.Width=roi(3);
            obj.Height=roi(4);
        end
    end
    
end
