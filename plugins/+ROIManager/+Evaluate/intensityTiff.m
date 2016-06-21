classdef intensityTiff<interfaces.SEEvaluationProcessor
    methods
        function obj=intensityTiff(varargin)        
                obj@interfaces.SEEvaluationProcessor(varargin{:});
        end
        function out=run(obj,p)
            % at some point: use shift_xy to fit with offset
            if ~isempty(p.layer2_)
            shiftx=p.layer2_.shiftxy_min;
            shifty=p.layer2_.shiftxy_max;
            else
                shiftx=0;
                shifty=0;
            end
            

            pos=obj.site.pos;
            pos(1)=pos(1)-shiftx;
            pos(2)=pos(2)-shifty;
            if p.posfromeval
                fieldx=p.dxfield;
                fieldy=p.dyfield;
                dx=eval(['obj.site.evaluation.' fieldx]);
                dy=eval(['obj.site.evaluation.' fieldy]);
                pos(1)=pos(1)+dx;
                pos(2)=pos(2)+dy;
               
            end
            
            file=obj.locData.files.file(obj.site.info.filenumber);
            pixcam=file.info.pixsize*1000;
            roi=file.tif(1).info.roi;

            pospixr=round(pos(1:2)/pixcam)-roi(1:2)';
            
            sfit=round((p.roisize-1)/2);
            
            coim=file.tif(1).image(pospixr(2)-sfit:pospixr(2)+sfit,pospixr(1)-sfit:pospixr(1)+sfit,:);
        
            pr=round(pos/pixcam);
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
                
                [fitp(3),fitp(4)]=restrictcoordiantes(fitp(3),fitp(4),pos(1),pos(2),p.mind);
%                 if fitp(3)<pos(1), fitp(3)=min(fitp(3),pos(1)-p.mind); else fitp(3)=max(fitp(3),pos(1)+p.mind); end
%                 if fitp(4)<pos(2), fitp(4)=min(fitp(4),pos(2)-p.mind); else fitp(4)=max(fitp(4),pos(2)+p.mind); end
                
                fitim=dgaussforfit(fitp,X,Y,fixp);
                
                startp2=[fitp pos(1) pos(2)];
                fixp2=fixp;
                fitp2=doublegaussfit2(coim,rangex,rangey,startp2,fixp2);
                
%                 fitimfree=dgaussforfit2(fitp2,X,Y,fixp2);
                ax=obj.setoutput('image');
                
                imagesc(rangex,rangey,coim,'Parent',ax)
                 axis(ax,'equal')
                ax.NextPlot='add';
                plot(pos(1),pos(2),'ko','Parent',ax)
                 plot(fitp2(7),fitp2(8),'kx','Parent',ax)
                d=sqrt((fitp(3)-pos(1))^2+(fitp(4)-pos(2)).^2);
                sx=(rangex(end)-rangex(1))/2;
                dv=(fitp(3:4)-pos(1:2));
                mv=max(abs(dv));
                if mv<=sx
                 plot(fitp(3),fitp(4),'kx','Parent',ax)
                else
                    dv=dv/mv*(sx+pixcam);
                    plot(pos(1)+dv(1),pos(2)+dv(2),'k*','Parent',ax)
                    ax.XLim=[rangex(1)-pixcam rangex(end)+pixcam];
                    ax.YLim=[rangey(1)-pixcam rangey(end)+pixcam];
                end
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
                out.xcent=fitp2(7);
                out.ycent=fitp2(8);
        end
        function pard=guidef(obj)
            pard=guidef;
        end
    end
end

function pard=guidef
pard.inputParameters={'numberOfLayers','sr_layerson','se_cellfov','se_sitefov','se_siteroi','layer2_'};
pard.t1.object=struct('Style','text','String','ROI (pixels)');
pard.t1.position=[1,1];
pard.t1.Width=2;

pard.roisize.object=struct('Style','edit','String','7');
pard.roisize.position=[1,4];

pard.t2.object=struct('Style','text','String','minimal distance 2nd Gauss (nm)');
pard.t2.position=[2,1];
pard.t2.Width=3;

pard.mind.object=struct('Style','edit','String','100');
pard.mind.position=[2,4];

pard.t3.object=struct('Style','text','String','sigma first Gauss (nm)');
pard.t3.position=[3,1];
pard.t3.Width=2;

pard.sigmaG.object=struct('Style','edit','String','150');
pard.sigmaG.position=[3,4];
pard.plugininfo.type='ROI_Evaluate';

pard.posfromeval.object=struct('Style','checkbox','String','use position from evaluation field: site.evaluation.');
pard.posfromeval.position=[5,1];
pard.posfromeval.Width=4;
pard.t4.object=struct('Style','text','String','dx0');
pard.t4.position=[6,1];
pard.t4.Width=.5;
pard.t5.object=struct('Style','text','String','dy0');
pard.t5.position=[7,1];
pard.t5.Width=.5;

pard.dxfield.object=struct('Style','edit','String','CME2DRing.imfit.x0');
pard.dxfield.position=[6,1.5];
pard.dxfield.Width=3.5;
pard.dyfield.object=struct('Style','edit','String','CME2DRing.imfit.y0');
pard.dyfield.position=[7,1.5];
pard.dyfield.Width=3.5;

% 
% pard.fitxy.object=struct('Style','checkbox','String','Fit also position of DL peaks');
% pard.fitxy.position=[6,1];
% pard.fitxy.Width=4;

pard.plugininfo.type='ROI_Evaluate';

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
% if x2<x1, x2=min(x2,x1-d); else x2=max(x2,x1+d); end
% if y2<y1, y2=min(y2,y1-d); else y2=max(y2,y1+d); end
[x2,y2]=restrictcoordiantes(x2,y2,x1,y1,d);
out=double(A1* exp( - ((X-x1).^2+(Y-y1).^2)/sigma1^2/2)+A2* exp( - ((X-x2).^2+(Y-y2).^2)/sigma2^2/2)+bg);

end



function fit=doublegaussfit2(img,rangex,rangey,startp,fixp)
[X,Y]=meshgrid(double(rangex),double(rangey));
lb=[0 0 -inf -inf 150 0 -inf -inf];ub=[inf inf inf inf 750 inf inf inf];
fit=lsqnonlin(@dgaussforfiterr2,double(startp),lb,ub,[],double(img),X,Y,fixp);

end
function out=dgaussforfiterr2(fitp,img,X,Y,fixp)
out=dgaussforfit2(fitp,X,Y,fixp)-img;
end
function out=dgaussforfit2(fitp,X,Y,fixp)
% fitp: A1 A2 x2 y2 sigma2 bg
%fixp x1, y1, sigma1

A1=fitp(1); A2=fitp(2); x2=fitp(3); y2=fitp(4); sigma2=fitp(5); bg=fitp(6);

d=fixp(4);
x1=fitp(7);y1=fitp(8); sigma1=fixp(3);
%adjsut x1, y1
% if x2<x1, x2=min(x2,x1-d); else x2=max(x2,x1+d); end
% if y2<y1, y2=min(y2,y1-d); else y2=max(y2,y1+d); end
[x2,y2]=restrictcoordiantes(x2,y2,x1,y1,d);
out=double(A1* exp( - ((X-x1).^2+(Y-y1).^2)/sigma1^2/2)+A2* exp( - ((X-x2).^2+(Y-y2).^2)/sigma2^2/2)+bg);

end

function [x2,y2]=restrictcoordiantes(x2,y2,x1,y1,d)
dvec=[x2-x1;y2-y1];
dist=norm(dvec);
% dist=sqrt((x2-x1)^2+(y2-y1)^2);
if dist<d
    xnew=[x1;x2]+dvec/dist*d;
    x2=xnew(1);
    y2=xnew(2);

end
end