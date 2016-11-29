classdef selectCoordinates2Color<interfaces.SEEvaluationProcessor
    properties
        hax1
        hax2
        hax3
        line
    end
    methods
        function obj=selectCoordinates2Color(varargin)        
                obj@interfaces.SEEvaluationProcessor(varargin{:});
        end
        function out=run(obj,p)
            obj.line=p.lineselect.selection;
            %  use shift_xy to fit with offset
            ax1=obj.setoutput('positions',true);
            hp=ax1.Parent;
%             delete(hp.Children);
%             ax1=axes(hp);
            ax2=axes(hp);
            ax3=axes(hp);
            subplot(1,3,1,ax1);
            subplot(1,3,2,ax2);
            subplot(1,3,3,ax3);
            site=obj.site;
            im1=site.image.layers(1).images.finalImages;
            if length(site.image.layers)<2
                ind=1;
            else
                ind=2;
            end
            im2=site.image.layers(ind).images.finalImages;
            
            range=[-p.se_sitefov p.se_sitefov]/2;
            imagesc(ax1,im1.rangex/1000,im1.rangey/1000,im1.image);
            axis(ax1,'equal');
            imagesc(ax2,im2.rangex/1000,im2.rangey/1000,im2.image);
            axis(ax2,'equal');
            imagesc(ax3,site.image.rangex,site.image.rangey,site.image.image);
            axis(ax3,'equal');
            pl=site.annotation.(obj.line).pos;
            if sum(pl)==0
                xp=mean(site.image.rangex);
                yp=mean(site.image.rangey);
                pl=[xp, yp;xp,yp];
                obj.site.annotation.(obj.line).pos=pl;
            end
            obj.hax1=impoint(ax1,pl(1,:));
            

            
            obj.hax1.addNewPositionCallback(@(p) obj.posconstraint(p,1));
            obj.hax2=impoint(ax2,pl(2,:));
            obj.hax2.addNewPositionCallback(@(p) obj.posconstraint(p,2));
            obj.hax3=imline(ax3,pl);
            obj.hax3.addNewPositionCallback(@(p) obj.posconstraint(p,3));
            out=[];

        end
        function pard=guidef(obj)
            pard=guidef;
        end
        function posconstraint(obj,pos,caller)
            obj.line
            linepos=obj.site.annotation.(obj.line).pos;
            
            switch caller
                case 1
                    linepos(1,:)=pos;
                case 2
                    linepos(2,:)=pos;
                case 3
                    linepos=pos;
            end
               
            obj.hax1.setPosition(linepos(1,:));
            obj.hax2.setPosition(linepos(2,:));
            obj.hax3.setPosition(linepos);
            obj.site.annotation.(obj.line).pos=linepos;
            
%             obj.site.annotation.line1.value=sqrt(sum((linepos(2,:)-linepos(1,:)).^2));
%             obj.site.annotation.line1.length=obj.site.annotation.line1.value;
        end
    end
end

function pard=guidef
pard.inputParameters={'numberOfLayers','sr_layerson','se_cellfov','se_sitefov','se_siteroi','layer2_'};

pard.lineselect.object=struct('Style','popupmenu','String',{{'line1','line2'}});
pard.lineselect.position=[1,1];
pard.lineselect.Width=4;

pard.plugininfo.type='ROI_Evaluate';

end


