classdef intensityTiff<interfaces.SEEvaluationProcessor
    methods
        function obj=intensityTiff(varargin)        
                obj@interfaces.SEEvaluationProcessor(varargin{:});
        end
        function out=run(obj,p)
            % at some point: use shift_xy to fit with offset
            shiftx=p.layer2_.shiftxy_min;
            shifty=p.layer2_.shiftxy_max;
            
            pos=obj.site.pos;
            pos(1)=pos(1)-shiftx;
            pos(2)=pos(2)-shifty;
            file=obj.locData.files.file(obj.site.info.filenumber);
            pixcam=file.info.pixsize*1000;
            roi=file.tif(1).info.roi;

            pospixr=floor(pos(1:2)/pixcam)-roi(1:2)';
            
            sfit=round((p.roisize-1)/2);
            
            coim=file.tif(1).image(pospixr(2)-sfit:pospixr(2)+sfit,pospixr(1)-sfit:pospixr(1)+sfit,:);
        
            pr=floor(pos/pixcam)+.5;
            rangex=(pr(1)-sfit:pr(1)+sfit)*pixcam;
            rangey=(pr(2)-sfit:pr(2)+sfit)*pixcam;
                
                
                fixp=[pos(1),pos(2),p.sigmaG,p.mind];
                [X,Y]=meshgrid(rangex,rangey);
                
                startim=dgaussforfit([max(coim(:)) 0 0 0 0 0],X,Y,[pos(1),pos(2), 250, p.mind]);
                [maxrem,ind]=max(coim(:)-startim(:));
                [my,mx]=ind2sub(size(startim),ind);
       
                startp=[max(coim(:)) maxrem rangex(mx) rangey(my) 150 0];
                startim=dgaussforfit(startp,X,Y, fixp);
          
                fitp=doublegaussfit(coim,rangex,rangey,startp,fixp);   
                
                
                if fitp(3)<pos(1), fitp(3)=min(fitp(3),pos(1)-p.mind); else fitp(3)=max(fitp(3),pos(1)+p.mind); end
                if fitp(4)<pos(2), fitp(4)=min(fitp(4),pos(2)-p.mind); else fitp(4)=max(fitp(4),pos(2)+p.mind); end
                
                fitim=dgaussforfit(fitp,X,Y,fixp);
                ax=obj.setoutput('image');
                
                imagesc(rangex,rangey,coim,'Parent',ax)
                ax.NextPlot='add';
                plot(pos(1),pos(2),'ko','Parent',ax)
                plot(fitp(3),fitp(4),'kx','Parent',ax)
                d=sqrt((fitp(3)-pos(1))^2+(fitp(4)-pos(2)).^2);
                ttxt=['A1=' num2str(fitp(1),'%5.0f') ', A2=' num2str(fitp(2),'%5.0f') ', d=' num2str(d,'%5.0f')];
                title(ttxt,'Parent',ax)
                axis(ax,'equal')
                ax.NextPlot='replace';
                ax2=obj.setoutput('fit');
                
                imagesc(rangex,rangey,fitim,'Parent',ax2)
                title(fitp(1:2),'Parent',ax2)
                axis(ax2,'equal')
                ax3=obj.setoutput('residuals');
                
                imagesc(rangex,rangey,coim-fitim,'Parent',ax3) 
                axis(ax3,'equal')
                
                ax4=obj.setoutput('startim');
                
                
                imagesc(rangex,rangey,startim,'Parent',ax4)   
                axis(ax4,'equal')
                out.Amplitude1=fitp(1);
                out.Amplitude2=fitp(2);
        end
        function pard=pardef(obj)
            pard=pardef;
        end
    end
end

function pard=pardef
pard.inputParameters={'numberOfLayers','sr_layerson','se_cellfov','se_sitefov','se_siteroi','layer2_'};
pard.t1.object=struct('Style','text','String','ROI (pixels)');
pard.t1.position=[1,1];
pard.t1.Width=2;

pard.roisize.object=struct('Style','edit','String','7');
pard.roisize.position=[1,3];

pard.t2.object=struct('Style','text','String','minimal distance 2nd Gauss (nm)');
pard.t2.position=[2,1];
pard.t2.Width=2;

pard.mind.object=struct('Style','edit','String','100');
pard.mind.position=[2,3];

pard.t3.object=struct('Style','text','String','sigma first Gauss (nm)');
pard.t3.position=[3,1];
pard.t3.Width=2;

pard.sigmaG.object=struct('Style','edit','String','150');
pard.sigmaG.position=[3,3];

end



function fit=doublegaussfit(img,rangex,rangey,startp,fixp)
[X,Y]=meshgrid(double(rangex),double(rangey));
lb=[0 0 -inf -inf 150 0];ub=[inf inf inf inf 750 inf];
fit=lsqnonlin(@dgaussforfiterr,double(startp),lb,ub,[],double(img),X,Y,fixp);

end
function out=dgaussforfiterr(fitp,img,X,Y,fixp)
out=dgaussforfit(fitp,X,Y,fixp)-img;
end
function out=dgaussforfit(fitp,X,Y,fixp)
% fitp: A1 A2 x2 y2 sigma2 bg
%fixp x1, y1, sigma1

A1=fitp(1); A2=fitp(2); x2=fitp(3); y2=fitp(4); sigma2=fitp(5); bg=fitp(6);

d=fixp(4);
x1=fixp(1);y1=fixp(2); sigma1=fixp(3);
%adjsut x1, y1
if x2<x1, x2=min(x2,x1-d); else x2=max(x2,x1+d); end
if y2<y1, y2=min(y2,y1-d); else y2=max(y2,y1+d); end
out=double(A1* exp( - ((X-x1).^2+(Y-y1).^2)/sigma1^2/2)+A2* exp( - ((X-x2).^2+(Y-y2).^2)/sigma2^2/2)+bg);

end