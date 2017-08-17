classdef BALM_fibril_growth<interfaces.SEEvaluationProcessor
    methods
        function obj=BALM_fibril_growth(varargin)        
                obj@interfaces.SEEvaluationProcessor(varargin{:});
        end
        function out=run(obj,p)
            out=[];
            lw=200;
            locs=obj.getLocs({'locprecnm','PSFxnm','xnm','ynm','frame'},'layer',1,'size',p.se_siteroi);  
            pol=obj.site.annotation.rotationpos.pos;
            
            dpol=pol(2,:)-pol(1,:);
            alpha=-atan2(dpol(1),dpol(2));

            len=sqrt(sum(dpol.^2))*1000/2;
            midp=mean(pol,1)*1000;
            [xr,yr]=rotcoord(locs.xnm-midp(1),locs.ynm-midp(2),alpha);

            indb=abs(xr)>lw|abs(yr)>len;
            indg=~indb;
            ynmline=xr(indg);xnmline=yr(indg);
            frame=locs.frame(indg);
            
            %segment filament, only keep locs in filament
            pixrec=10; %nm
            xrange=-len-pixrec:pixrec:len+pixrec;
            yrange=-lw-pixrec:pixrec:lw+pixrec;
            him=histcounts2(xnmline, ynmline,xrange,yrange);
            sigma=0.5;
            h=fspecial('gaussian',12,sigma);
            max1=mean(him(:));
            
            himf=fibermetric(him,15,'StructureSensitivity',max1/3);
            himf=filter2(h,himf);
             max1=mean(himf(:));
%             max1=1/sigma^2/pi/2;
%             max1=4*max(h(:));
            
            imbw=himf>max1;
%             figure(89);subplot(2,2,1);imagesc(himf);
            imbw=bwareaopen(imbw,round(sum(imbw(:))/4));
%              subplot(2,2,2);imagesc(imbw);
            
             seb=strel('disk',3,4);
%              ses=strel('disk',1,4);
%              imbw=imerode(imbw,ses);
%              imbw=bwareaopen(imbw,round(sum(imbw(:))/4));
            imbw=imdilate(imbw,seb);
            imbw=imerode(imbw,seb);
             h=obj.setoutput('images');
             himp=him/max(him(:));
             imagesc(h,horzcat(himp,himf,double(imbw)))
%             
%            imbw=imdilate(imbw,ses);
%             imbw=imerode(imbw,ses);
%             imbw=imopen(imbw,seb);
%              imbw=imerode(imbw,ses);
%               imbw=imerode(imbw,ses);

           
           
         
%             imbw=imsegfmm(himf,imbw,0.01);
%             imbw=imdilate(imbw,ones(3));
%             bwac=activecontour(himf,true(size(himf)));
%             subplot(2,2,3);imagesc(imbw)
            xr=round((xnmline+len)/pixrec)+1;
            yr=round((ynmline+lw)/pixrec)+1;
            indlin=sub2ind(size(imbw),xr,yr);
            indgood=imbw(indlin);
            
            %
%             h2=fspecial('gaussian',12,1.5);
            
 
            xr=-len:p.dx:len;
            fr=1:p.df:max(frame(indgood));
            him=histcounts2(frame(indgood),xnmline(indgood),fr,xr);
%             hxf=filter2(h2,him);
%             mf=max(h2(:));
%             hxfb=hxf>mf*1.1;
            
            
            co=quantile(him(:),0.999);
            him(him>co)=co;
            h=obj.setoutput('kimograph');
            imagesc(h,xr,fr,(him))
            xlabel(h,'xnm')
            ylabel(h,'frame')
            



        end
        function pard=guidef(obj)
            pard=guidef;
        end
    end
end

function pard=guidef
pard.dxt.object=struct('Style','text','String','dx (nm), dt (frames)');
pard.dxt.position=[1,1];
pard.dxt.Width=2;
pard.dx.object=struct('Style','edit','String','10');
pard.dx.position=[2,1];
pard.df.object=struct('Style','edit','String','1000');
pard.df.position=[2,2];

% pard.dxt.Width=3;
pard.inputParameters={'numberOfLayers','sr_layerson','se_cellfov','se_sitefov','se_siteroi'};
pard.plugininfo.type='ROI_Evaluate';
end