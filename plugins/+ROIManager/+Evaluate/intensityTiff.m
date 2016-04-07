classdef intensityTiff<interfaces.SEEvaluationProcessor
    methods
        function obj=intensityTiff(varargin)        
                obj@interfaces.SEEvaluationProcessor(varargin{:});
        end
        function out=run(obj,p)
            pos=obj.site.pos;
            file=obj.locData.files.file(obj.site.info.filenumber);
            pixcam=file.info.pixsize*1000;
            roi=file.tif(1).info.roi;

            pospixr=floor(pos(1:2)/pixcam)-roi(1:2)';
            
            sfit=(7-1)/2;
            
            coim=file.tif(1).image(pospixr(2)-sfit:pospixr(2)+sfit,pospixr(1)-sfit:pospixr(1)+sfit,:);
        
            pr=floor(pos/pixcam)+.5;
            rangex=(pr(1)-sfit:pr(1)+sfit)*pixcam;
            rangey=(pr(2)-sfit:pr(2)+sfit)*pixcam;
%                 figure(88);imagesc(rangex,rangey,coim)
                
                
                fixp=[pos(1),pos(2),150];
                [X,Y]=meshgrid(rangex,rangey);
                
                startim=dgaussforfit([max(coim(:)) 0 0 0 0 0],X,Y,[pos(1),pos(2) 250]);
                [maxrem,ind]=max(coim(:)-startim(:));
                [my,mx]=ind2sub(size(startim),ind);
       
                startp=[max(coim(:)) maxrem rangex(mx) rangey(my) 250 0];
                startim=dgaussforfit(startp,X,Y, fixp);
%                 figure(89); imagesc(rangex,rangey,dgaussforfit(startp,X,Y,fixp));
          
                fitp=doublegaussfit(coim,rangex,rangey,startp,fixp);      
%                 fitp(3:4)-fixp(1:2)
                fitim=dgaussforfit(fitp,X,Y,fixp);
                ax=obj.setoutput('image');
                axis(ax,'equal')
                imagesc(rangex,rangey,coim,'Parent',ax)
                ax2=obj.setoutput('fit');
                axis(ax2,'equal')
                imagesc(rangex,rangey,fitim,'Parent',ax2)
                title(fitp(1:2),'Parent',ax2)
                ax3=obj.setoutput('residuals');
                axis(ax3,'equal')
                imagesc(rangex,rangey,coim-fitim,'Parent',ax3)      
                ax4=obj.setoutput('startim');
                axis(ax4,'equal')
                imagesc(rangex,rangey,startim,'Parent',ax4)      
                out.Amplitude1=fitp(1);
                out.Amplitude2=fitp(2);
        end
        function pard=pardef(obj)
            pard=pardef;
        end
    end
end

function pard=pardef
pard.inputParameters={'numberOfLayers','sr_layerson','se_cellfov','se_sitefov','se_siteroi'};
end



function fit=doublegaussfit(img,rangex,rangey,startp,fixp)
[X,Y]=meshgrid(rangex,rangey);
lb=[0 0 -inf -inf 150 0];ub=[inf inf inf inf 750 inf];
fit=lsqnonlin(@dgaussforfiterr,startp,lb,ub,[],img,X,Y,fixp);

end
function out=dgaussforfiterr(fitp,img,X,Y,fixp)
out=dgaussforfit(fitp,X,Y,fixp)-img;
end
function out=dgaussforfit(fitp,X,Y,fixp)
% fitp: A1 A2 x2 y2 sigma2 bg
%fixp x1, y1, sigma1

A1=fitp(1); A2=fitp(2); x2=fitp(3); y2=fitp(4); sigma2=fitp(5); bg=fitp(6);

d=100;
x1=fixp(1);y1=fixp(2); sigma1=fixp(3);
%adjsut x1, y1
if x2<x1, x2=min(x2,x1-d); else x2=max(x2,x1+d); end
if y2<y1, y2=min(y2,y1-d); else y2=max(y2,y1+d); end
out=A1* exp( - ((X-x1).^2+(Y-y1).^2)/sigma1^2/2)+A2* exp( - ((X-x2).^2+(Y-y2).^2)/sigma2^2/2)+bg;

end