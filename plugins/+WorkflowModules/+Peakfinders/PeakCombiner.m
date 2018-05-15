classdef PeakCombiner<interfaces.WorkflowModule
    properties (Access=private)
        transform
    end
    methods
        function obj=PeakCombiner(varargin)
            obj@interfaces.WorkflowModule(varargin{:});
             obj.inputParameters={'loc_loc_filter_sigma','EMon','loc_ROIsize','loc_fileinfo'};
        end
        function pard=guidef(obj)
            pard=guidef(obj);
        end
        function initGui(obj)
            initGui@interfaces.WorkflowModule(obj);
%             obj.guihandles.cutoffmode.Callback={@cutoffmode_callback,obj};
%             obj.guihandles.cutoffvalue.Callback={@cutoffvalue_callback,obj};
%              obj.guihandles.peakfindmethod.Callback={@peakfindmethod_callback,obj};
%              cutoffvalue_callback(0,0,obj)
        end
        function prerun(obj,p)
            l=load(p.Tfile);
            if isfield(l,'transformation')
                obj.transform=l.transformation;  
            else
                obj.transform=l.SXY(1).cspline.global.transformation;
            end
            obj.setPar('loc_globaltransform',obj.transform);
%             cutoffvalue_callback(0,0,obj)
        end
        function dato=run(obj,data,p)
            if isempty(data.data)
                dato=data;
                return;
            end
            roi=p.loc_fileinfo.roi;
%             roi=zeros(1,4);
            maxima=data.data;%get;
            transform=obj.transform;
            if isfield(transform.tinfo,'cam_pixelsize_nm')
                pixelsize=transform.tinfo.cam_pixelsize_nm;
            else
                pixelsize=p.loc_fileinfo.cam_pixelsize_um*1000;
            end
            if isfield(transform.tinfo,'units')&&strcmp(transform.tinfo.units,'pixels')
                pixelsize=[1 1];
            end
            xnm=(maxima.x+roi(1))*pixelsize(1);
            ynm=(maxima.y+roi(2))*pixelsize(end);
            indref=obj.transform.getRef(xnm,ynm);
            xnmref=xnm(indref);
            ynmref=ynm(indref);
            Nref=maxima.intensity(indref).^2;
            Nt=maxima.intensity(~indref).^2;
            [xt,yt]=transform.transformCoordinatesInv(xnm(~indref),ynm(~indref));
            [iA,iB,uiA,uiB]=matchlocs(xnmref,ynmref,xt,yt,[0 0],pixelsize(1)*4);
            try
            xallr=vertcat(xnmref(uiA), xt(uiB'),((xnmref(iA).*Nref(iA)+xt(iB).*Nt(iB))./(Nref(iA)+Nt(iB))));
            yallr=vertcat(ynmref(uiA), yt(uiB'),((ynmref(iA).*Nref(iA)+yt(iB).*Nt(iB))./(Nref(iA)+Nt(iB))));
            catch err
                disp('error')
            end
%             
%             xallr=vertcat(xnmref(uiA), xnmref(iA));
%             yallr=vertcat(ynmref(uiA), ynmref(iA));
            
            
            xallrpix=round(xallr/pixelsize(1));
            yallrpix=round(yallr/pixelsize(end));
            
            %for challenge: allow for localizations close to rim. Correct
            %both reference and target if too far outside
            %now implemented only for right-left
            if contains(p.Tfile,'Challenge')
                dn=(p.loc_ROIsize-1)/2;
                xallrpix=min(max(xallrpix,dn+1),roi(3)/2-dn);
                yallrpix=min(max(yallrpix,dn+1),roi(4)-dn);
            end

             
            [xallt,yallt]=transform.transformCoordinatesFwd(xallrpix*pixelsize(1),yallrpix*pixelsize(end));
            
            xalltpix=round(xallt/pixelsize(1));
            yalltpix=round(yallt/pixelsize(end)); 
            

            
            indgood=xalltpix>roi(1)&yalltpix>roi(2)&xalltpix<=roi(1)+roi(3)&yalltpix<=roi(2)+roi(4);
            
            dx=xallt(indgood)/pixelsize(1)-xalltpix(indgood);
            dy=yallt(indgood)/pixelsize(end)-yalltpix(indgood);
            
            xallrpix=xallrpix(indgood)-roi(1);
            yallrpix=yallrpix(indgood)-roi(2);
            xalltpix=xalltpix(indgood)-roi(1);
            yalltpix=yalltpix(indgood)-roi(2);
            
        
            xout=vertcat(xallrpix',xalltpix');
            yout=vertcat(yallrpix',yalltpix');
            dxout=vertcat(0*dx',dx');
            dyout=vertcat(0*dy',dy');
            idout=vertcat(1:sum(indgood),1:sum(indgood));
            
            maxout.x=xout(:);
            maxout.y=yout(:);
            maxout.xfound=maxima.x;
            maxout.yfound=maxima.y;
            maxout.intensityfound=maxima.intensity;
            maxout.ID=idout(:);
            maxout.dx=dxout(:);
            maxout.dy=dyout(:);
            dato=data;
            dato.data=maxout;
        end
        function loadbutton(obj,a,b)
            fn=obj.guihandles.Tfile.String;
            path=fileparts(fn);
            if ~exist(path,'file')
                fn=[fileparts(obj.getPar('loc_outputfilename')) filesep '*.mat'];
            end
            [f,path]=uigetfile(fn,'Select transformation file _T.mat');
            if f
                Tload=load([path f],'transformation');
                if ~isfield(Tload,'transformation')
                    msgbox('could not find transformation in file. Load other file?')
                end
                
                obj.guihandles.Tfile.String=[path f];
                obj.setPar('transformationfile',[path f]);
            end      
        end  
    end
end



function pard=guidef(obj)
pard.Tfile.object=struct('Style','edit','String','settings/temp/temp_T.mat');
pard.Tfile.position=[1,1];
pard.Tfile.Width=3;

pard.loadbutton.object=struct('Style','pushbutton','String','load T','Callback',@obj.loadbutton);
pard.loadbutton.position=[1,4];

pard.syncParameters={{'transformationfile','Tfile',{'String'}}};
% pard.cutoffmode.object=struct('Style','popupmenu','String',{{'dynamic (factor)','probability (p<1)','absolute (photons)'}},'Value',2);
% pard.cutoffmode.position=[1,1];
% pard.cutoffmode.Width=1.5;
% pard.cutoffmode.TooltipString=sprintf('How to determine the cutoff: \n Dynamic: use the distribution of pixel intensity to estimate likely localizations. Factor: adjust sensitivity. \n Probability: use probabilistic model (SimpleSTORM) to determine the likelyhood for pixel being localization. \n Absolut: Pixel intensity in normalized image. \n Choose display=normalized to read out thes normalized values.');
% pard.cutoffmode.Optional=true;

pard.plugininfo.type='WorkflowModule'; 
pard.plugininfo.description='Performs peak-finding according either to local maximum finding  or to non-maximum suppression (A. Neubeck and L. Van Gool, ?Efficient non-maximum suppression,? presented at the 18th International Conference on Pattern Recognition, Vol 3, Proceedings, 10662 LOS VAQUEROS CIRCLE, PO BOX 3014, LOS ALAMITOS, CA 90720-1264 USA, 2006, pp. 850?855.).';
end