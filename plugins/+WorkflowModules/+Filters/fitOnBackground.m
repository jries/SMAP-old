classdef fitOnBackground<interfaces.WorkflowModule
    properties
        fitonbg

    end
    methods
        function obj=fitOnBackground(varargin)
            obj@interfaces.WorkflowModule(varargin{:})
            obj.inputChannels=2; 
        end
        function pard=pardef(obj)
            pard=pardef;
        end
        function initGui(obj)
            initGui@interfaces.WorkflowModule(obj);
            obj.setInputChannels(2,'frame');
        end
        function prerun(obj,p)
            obj.fitonbg=p.loc_fitonbg;
%              obj.preview=obj.getPar('loc_preview');
        end
        function output=run(obj,data,p)
            output=[];
       
            if  ~isempty(data{1}.data)
                
                image=data{1}.data.img;%get; 
                bg=data{2}.data.img;%get;
                if obj.fitonbg
                    s=size(image);
                    imbgc=zeros(s,'like',image);
                    bgbgc=zeros(s,'like',image);
                    if length(s)==2
                        s(3)=1;
                    end

                    for k=1:s(3)
                        bgh=bg(:,:,k);
                        imh=image(:,:,k);
                        indnan=isnan(bgh);
                        meanbg=mean(bgh(~indnan));
                        imbgc(:,:,k)=imh-bgh+meanbg;
                        bgbgc(:,:,k)=meanbg;
                    end
                else
                    imbgc=image;
                    bgbgc=bg;
                end

                dato=data{1};%{1}.copy;
                dato.data.img=imbgc;%set(imnorm);
                obj.output(dato,1)

                dato2=data{2};%{1}.copy;
                dato2.data.img=bgbgc;%set(imnorm);
                obj.output(dato2,2)
            else 
                obj.output(data{1},1);
                obj.output(data{2},2);
            end
        end
        
    end
end


function pard=pardef
pard.loc_fitonbg.object=struct('Style','checkbox','String','fit on BG');
pard.loc_fitonbg.position=[1,1];
pard.loc_fitonbg.Width=1;
end