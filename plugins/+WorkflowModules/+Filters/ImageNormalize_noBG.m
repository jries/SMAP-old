classdef ImageNormalize_noBG<interfaces.WorkflowModule
    properties
        preview

    end
    methods
        function obj=ImageNormalize_noBG(varargin)
            obj@interfaces.WorkflowModule(varargin{:})
            obj.inputChannels=1; 
        end
        function pard=guidef(obj)
            pard.plugininfo.type='WorkflowModule';
            pard.plugininfo.description='Converts photons into a probability map. According to: [1]	U. Koethe, F. Herrmannsdoerfer, I. Kats, and F. A. Hamprecht, SimpleSTORM: a fast, self-calibrating reconstruction algorithm for localization microscopy,HISTOCHEMISTRY AND CELL BIOLOGY, pp. 1-15, Apr. 2014.';
        end
        function initGui(obj)
            initGui@interfaces.WorkflowModule(obj);
        end
        function prerun(obj,p)          
            obj.preview=obj.getPar('loc_preview');          
        end
        function dato=run(obj,data,p)
            
            if  ~isempty(data.data)
                image=data.data;%get;  
                imnorm=poissonNormalize(image);
                dato=data;%{1}.copy;
                dato.data=imnorm;%set(imnorm);
%                 if obj.preview
%                     if data.frame==obj.getPar('loc_previewframe')
%                         drawimage(obj,imnorm,image,0)
%                     else
%                         dato=[];
%                     end
%                 end

            else 
                dato=data;
            end     
        end
        

    end
end


% 
% function drawimage(obj,imnorm,img,bg)
% outputfig=obj.getPar('loc_outputfig');
% if ~isvalid(outputfig)
%     outputfig=figure(209);
%     obj.setPar('loc_outputfig',outputfig);
% end
% 
% outputfig.Visible='on';
% draw=true;
% switch obj.getPar('loc_previewmode').Value
%     case 1 %image-bg
%         imd=img-bg;
%     case 2%image
%         imd=img;
% %     case 3 %norm
% %         imd=imnorm;
%     case 4 %bg
%         imd=bg;
%     otherwise 
%         draw=false;
% end
%         
% if draw
% figure(outputfig)
% hold off
% imagesc(imd);
% colormap jet
% colorbar;
% axis equal
% end
% end

function out=poissonNormalize(in)
% out=real(2*sqrt(in+3/8));
in(in<-0.3750)=0;
out=(2*sqrt(in+0.3750));
end
