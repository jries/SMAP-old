classdef Viewer3DV01<interfaces.DialogProcessor
    properties
        axis
        timer
        currentimage
        commandfig
        recpar={};
        stereopar
    end
    methods
        function obj=Viewer3DV01(varargin)        
            obj@interfaces.DialogProcessor(varargin{:}) ;
            obj.inputParameters=[renderSMAP drawerSMAP displayerSMAP];
            obj.inputParameters{end+1}='sr_roihandle';
            obj.inputParameters{end+1}='linewidth_roi';
            obj.inputParameters{end+1}='layers';
            obj.inputParameters{end+1}='numberOfLayers';
            obj.inputParameters=unique(obj.inputParameters);
             obj.showresults=false;
        end
        function makeGui(obj)
            makeGui@interfaces.DialogProcessor(obj);
            h=obj.guihandles;           
        end
        
        function showpanel_callback(obj,a,b)
            if isempty(obj.commandfig)||~isvalid(obj.commandfig)
                
                fp=getParentFigure(obj.handle);
                posfig=fp.Position;
                posfig(1)=posfig(1)+posfig(3);
                posfig(3:4)=[200,200];
                obj.commandfig=figure('MenuBar','none','Toolbar','none','Position',posfig);
            end
            figure(obj.commandfig);
             h.ptranslation=makepanel([0 0.5 0.5 0.5],'translation','');
            h.protation=makepanel([0.5 0.5 0.5 0.5],'rot (command/strg)','command');
            h.pzoom=makepanel([0.5 0.0 0.5 0.5],'zoom (alt)','alt');

            function h=makepanel(pos,title,modifier)
                h=uipanel('Parent',obj.commandfig,'Units','normalized','Position',pos,'Title',title);
                h.Units='normalized';
                uicontrol('Parent',h,'Units','normalized','Position',[1 1 1 1]/3,'String','0','Callback',{@obj.keypress,struct('Modifier',modifier,'Key','0','Character','0')})
                uicontrol('Parent',h,'Units','normalized','Position',[0 1 1 1]/3,'String','<-','Callback',{@obj.keypress,struct('Modifier',modifier,'Key','leftarrow','Character','')})
                uicontrol('Parent',h,'Units','normalized','Position',[2 1 1 1]/3,'String','->','Callback',{@obj.keypress,struct('Modifier',modifier,'Key','rightarrow','Character','')})
                uicontrol('Parent',h,'Units','normalized','Position',[1 2 1 1]/3,'String','^','Callback',{@obj.keypress,struct('Modifier',modifier,'Key','uparrow','Character','')})
                uicontrol('Parent',h,'Units','normalized','Position',[1 0 1 1]/3,'String','v','Callback',{@obj.keypress,struct('Modifier',modifier,'Key','downarrow','Character','')})   
                uicontrol('Parent',h,'Units','normalized','Position',[0 2 1 1]/3,'String','x (,)','Callback',{@obj.keypress,struct('Modifier',modifier,'Key','comma','Character',',')})
                uicontrol('Parent',h,'Units','normalized','Position',[2 2 1 1]/3,'String','(.)','Callback',{@obj.keypress,struct('Modifier',modifier,'Key','period','Character','.')})                  
            end
        end
        function initGui(obj)

        end
        function out=run(obj,p)
            out=[];
            obj.addSynchronization('sr_roiposition',[],[],{@obj.redraw});
            
            if isempty(obj.axis)||~isvalid(obj.axis)
                figure;
                obj.axis=gca;
            end
            obj.axis.Units='normalized';
            obj.axis.Position=[0.05 0.05 .95 .9];
            
             set(obj.axis,'NextPlot','replacechildren','PickableParts','all','Units','pixels')
            fig=getParentFigure(obj.axis);
             axis(obj.axis,'tight');
             axis(obj.axis,'equal');
             axis(obj.axis,'ij');
            
             set(fig,'WindowKeyPressFcn',{@obj.keypress,[]})
             obj.timer=uint64(0);
             obj.redraw
        end

        function pard=guidef(obj)
            pard=guidef(obj);
        end
        function keypress(obj,a,d2,data)
            if isempty(data)
                data=d2;
            end
           %1.up, 2.down, 3.left, 4.right, 5.back, 6.front, 0.reset
            switch data.Character
                case {'w','8',30}
                    dir=1;
                case {'s','2',31}
                    dir=2;
                case {'a','4',28}
                    dir=3;
                case {'d','6',29}
                    dir=4;
                case {',','q','7'}
                    dir=5;
                case {'.','e','9'}
                    dir=6;
                case {'0'}
                    dir=0;   
               
                otherwise 
                    return
            end
            
            p=obj.getGuiParameters;
            roih=obj.getPar('sr_roihandle');
            pos=roih.getPosition;
            roivec=pos(2,:)-pos(1,:);
            roivecp(2)=roivec(1);
            roivecp(1)=-roivec(2);
            step=0.1;
            stepl=0.3;
            dphi=pi/16;
            dtheta=pi/16;
            if any(strcmp(data.Modifier,'shift'))
                stepfac=0.2;
            else
                stepfac=1;
            end
%             data.Modifier
            if any(strcmp(data.Modifier,'command'))||any(strcmp(data.Modifier,'control'))
                %rotate
                phi=0;
                theta=p.theta;
                switch dir
                    
                    case 1
                       %tilt up down
                       theta=theta+dtheta*stepfac;
                       if theta>pi
                           theta=theta-2*pi;
                       end
                       if theta<-pi
                           theta=theta+2*pi;
                       end
                       phi=0;
                    case 2
                        phi=0;
                       theta=theta-dtheta*stepfac;
                       if theta>pi
                           theta=theta-2*pi;
                       end
                       
                       if theta<-pi
                           theta=theta+2*pi;
                       end
                    case 3
                       phi=dphi*stepfac;
                    case 4
                        phi=-dphi*stepfac;
                    case 0
                        if strcmp(data.Character,'0')
                        theta=0;
                        end
                end
                 mpos=mean(pos,1);
                [dx,dy]=rotcoord(roivec(1)/2,roivec(2)/2,phi);
                pos(1,1)=mpos(1)-dx;
                pos(2,1)=mpos(1)+dx;
                pos(1,2)=mpos(2)-dy;
                pos(2,2)=mpos(2)+dy;
                obj.setGuiParameters(struct('theta',theta));
            elseif any(strcmp(data.Modifier,'alt'))
                %change size
                switch dir
                    case 6
                        lw=obj.getPar('linewidth_roi');
                        lw2=lw*(1+step*stepfac);
                        obj.setPar('linewidth_roi',lw2);
                       %tilt up down
                    case 5
                        lw=obj.getPar('linewidth_roi');
                        lw2=lw*(1-step*stepfac);
                        obj.setPar('linewidth_roi',lw2);
                      
                    case 3
                        
                       pos(1,:)=pos(1,:)+roivec/2*step*stepfac;
                       pos(2,:)=pos(2,:)-roivec/2*step*stepfac;
                    case 4
                        pos(1,:)=pos(1,:)-roivec/2*step*stepfac;
                       pos(2,:)=pos(2,:)+roivec/2*step*stepfac;
                    case 1
                        po.zmin=p.zmin-stepfac*step*(p.zmax-p.zmin);
                        po.zmax=p.zmax+stepfac*step*(p.zmax-p.zmin);
                        obj.setGuiParameters(po);
                    case 2
                        po.zmin=p.zmin+stepfac*step*(p.zmax-p.zmin);
                        po.zmax=p.zmax-stepfac*step*(p.zmax-p.zmin);
                        obj.setGuiParameters(po);
                end
            else
                switch dir
                    case 6
                        lw=obj.getPar('linewidth_roi')/1000;
                        pos(1,:)=pos(1,:)+stepl*roivecp./norm(roivecp)*stepfac*lw;
                        pos(2,:)=pos(2,:)+stepl*roivecp./norm(roivecp)*stepfac*lw;
                    case 5
                        lw=obj.getPar('linewidth_roi')/1000;
                        pos(1,:)=pos(1,:)-stepl*roivecp./norm(roivecp)*stepfac*lw;
                        pos(2,:)=pos(2,:)-stepl*roivecp./norm(roivecp)*stepfac*lw;
                    case 3
                        pos(1,:)=pos(1,:)+step*roivec*stepfac;
                        pos(2,:)=pos(2,:)+step*roivec*stepfac;
                    case 4
                        pos(1,:)=pos(1,:)-step*roivec*stepfac;
                        pos(2,:)=pos(2,:)-step*roivec*stepfac;
                    case 1
                        po.zmin=p.zmin+stepfac*step*(p.zmax-p.zmin);
                        po.zmax=p.zmax+stepfac*step*(p.zmax-p.zmin);
                        obj.setGuiParameters(po);
                    case 2
                        po.zmin=p.zmin-stepfac*step*(p.zmax-p.zmin);
                        po.zmax=p.zmax-stepfac*step*(p.zmax-p.zmin);
                        obj.setGuiParameters(po);
                    case 0
                        if strcmp(data.Character,'0')
                        po.zmin=p.zmin-(p.zmax+p.zmin)/2;
                        po.zmax=p.zmax-(p.zmax+p.zmin)/2;
                        obj.setGuiParameters(po);
                        end
                end
            end
            roih.setPosition(pos);
        end
        
        
        function redraw(obj)
            
            p=obj.getAllParameters;
            stereo=p.stereo.Value>1;
            locCopy=obj.locData; %maybe not needed
            lo=logical(obj.getPar('sr_layerson'));
            layerson=find(lo);
            if sum(lo)==0
                return
            end
            nl=obj.getPar('numberOfLayers');
            g=locCopy.isgrouped(1:nl);
            
            gt=g(lo(1:nl));
            group=zeros(2,1);
            if any(gt==1)
                group(2)=1;
            end
            if any(gt==0)
                group(1)=1;    
            end
            
            for k=nl:-1:1
                renderfield{k}=p.(['layer' num2str(k) '_']).renderfield.selection;
            end
            
            zmean=(p.zmin+p.zmax)/2;
            if p.setpixelsize
                ph.sr_pixrec=p.pixrecset;
            else
                ph.sr_pixrec=p.sr_pixrec;
            end  
            

            ph.rangey=[p.zmin p.zmax];
            roih=obj.getPar('sr_roihandle');
            rpos=roih.getPosition;

            lr=sqrt(sum((rpos(2,:)-rpos(1,:)).^2));
            rx=[-lr/2 lr/2]*1000;
            ax=obj.axis;
            obj.timer=tic;
           
            ph.sr_roihandle=obj.getPar('sr_roihandle');       
            ph.rangex=rx;
            
            if group(1)
                [loc,indu,sortind]=getlocrot('ungrouped','inlayeru');                
            end
            
            if group(2)
                [locg,indg,sortindg]=getlocrot('grouped','inlayerg');
            end
           
            %transparency
            transparency.parameter=p.transparencypar;
            transparency.mode=p.transparencymode.Value;
            switch p.transparencymode.Value
                case 1 %MIP
                case 2 %transparent
                    transparency.parameter=p.transparencypar/ph.sr_pixrec(1)^2*10;
                case 3 %balls
            end
%             if p.settransparency
%             transparency=p.transparency/ph.sr_pixrec(1)^2*10;
%             else
%                 transparency=[];
%             end
            
            
            ph.sr_axes=[];
            for k=1:p.numberOfLayers
                pl=p.(['layer' num2str(k) '_']);
                if pl.layercheck
                    if length(obj.recpar)>=k
                        rp=obj.recpar{k};
                    else
                        rp=[];
                    end
                     pr=copyfields(copyfields(copyfields(p,pl),ph),rp);
                     if stereo
                         pr=getstereosettings(pr,1);
                         layer1(k).images=renderplotlayer(pr,1);
                         %same intensity scaling
                         pr.imax=layer1(k).images.finalImages.imax;
                         obj.currentimage.imax(k)=pr.imax;
                         pr.imaxtoggle=0;
                         pr=getstereosettings(pr,2);
                         layer2(k).images=renderplotlayer(pr,2);
                     else
                        layer(k).images=renderplotlayer(pr,0);
                        obj.currentimage.imax(k)=layer(k).images.finalImages.imax;
                     end
                end
            end
            if stereo
                srim1=displayerSMAP(layer1,pr);
                srim2=displayerSMAP(layer2,pr);
                switch p.stereo.Value
                    case 2
                        srim=srim1;
                        srim.image=srim1.image+srim2.image;
                    case 3
                        srim=assembleSideviews(srim2,srim1,p);
                    case 4
                        srim=assembleSideviews(srim1,srim2,p);
                end
            else
                srim=displayerSMAP(layer,pr);
            end
            obj.currentimage=copyfields(obj.currentimage,srim);
            imagesc(srim.rangex*1000,srim.rangey*1000,srim.image,'Parent',ax);
           drawnow limitrate 
           
           
           
            function srim=assembleSideviews(srim1,srim2,p)
                srim=srim1;
                mpar=p.stereomagnification;
                widths=obj.axis.Parent.Position(3);
                
                ar=obj.axis.Parent.Position(3)/obj.axis.Parent.Position(4);
                sx=ceil(round(widths*mpar/2))*2;
                sy=ceil(round(widths*mpar/ar/2))*2;
                image=zeros(sy,sx,3);
                s=size(srim1.composite);
                mpx=ceil(sx/2);
                mpy=ceil(sy/2);
                shalfimy=floor(s(1)/2);
                shalfimx=floor(s(2)/2);
                eyepx=round(obj.stereopar.eyemm/p.monitorpmm*mpar);
                mpxh=mpx-eyepx;
                shalfy=min(shalfimy,mpy);
                shalfx=min(shalfimx,eyepx);
%                 im1=srim1.composite;
                rim1=shalfimy-shalfy+1:shalfimy+shalfy;
                rim2=shalfimx+1-shalfx:shalfimx+shalfx;
                imcut1=srim1.composite(rim1,rim2,:);
                imcut2=srim2.composite(rim1,rim2,:);
                
                image(mpy-shalfy+1:mpy+shalfy,mpxh-shalfx+1:mpxh+shalfx,:)=imcut1;
                mpxh=mpx+eyepx;
                image(mpy-shalfy+1:mpy+shalfy,mpxh-shalfx+1:mpxh+shalfx,:)=imcut2;
%                 image(mpy-shalf(1)+1:mpy-shalf(1)+s(1),mpxh-shalf(2)+1:mpxh-shalf(2)+s(2),:)=imcut2;                
                
                fx=sx/s(2);fy=sy/s(1);
                srim.image=image;
                srim.rangex=srim1.rangex*fx;
                srim.rangey=srim1.rangey*fy;
            end
            function images=renderplotlayer(pr,stereochannel)
                if stereochannel>0
                    if pr.groupcheck
                        locg.x=locg.(['x' num2str(stereochannel)]);
                    else
                        loc.x=loc.(['x' num2str(stereochannel)]);
                    end
                end
                 if pr.groupcheck
                        indroi=locg.inlayerg{layerson(k)};
                        indh=(indroi(indg));
                        images.srimage=renderSMAP(locg,pr,k,indh(sortindg),transparency);
                 else 
                     indroi=loc.inlayeru{layerson(k)};
                     indh=(indroi(indu));
                     images.srimage=renderSMAP(loc,pr,k,indh(sortind),transparency);
                 end
                images.finalImages=drawerSMAP(images.srimage,pr);        
                
            end
            function [loc,indu,sortind]=getlocrot(grouping,inlayer)
                [loc,indu]=locCopy.getloc({'xnmline','ynmline','znm','locprecnm','locprecznm',renderfield{:},inlayer},'position','roi','grouping',grouping,'layer',layerson);             
                [yrot,depth]=rotcoord(loc.znm-zmean,loc.ynmline,p.theta);
                [sortdepth,sortind]=sort(depth);
%                 sortdepth=sortdepth-max(sortdepth); %reference point on plane
                if stereo
                    eyemm=35;
                    dplanemm=p.dplanemm;
                    pixmm=p.monitorpmm;
                    widthpix=obj.axis.Position(3);
%                     pixnm=ph.sr_pixrec;
                    widthnm=ph.rangex(2)-ph.rangex(1);
                    heightpix=obj.axis.Position(4);
                    widthmm=widthpix*pixmm;
                    eyenm=eyemm*widthnm/widthmm;
                    dplanenm=dplanemm*widthnm/widthmm;
                    xe1=eyenm;xe2=-eyenm;
                    
                    x=loc.xnmline(sortind);
                    loc.x1=(x-xe1)./(1-sortdepth/dplanenm)+xe1;
                    loc.x2=(x-xe2)./(1-sortdepth/dplanenm)+xe2;
                    
                    obj.stereopar=struct('eyemm',eyemm,'eyenm',eyenm,'widthnm',widthnm,'widthmm',widthmm,'widthpix',widthpix,'heightpix',heightpix);
                else
                    loc.x=loc.xnmline(sortind);
                end
    
                
                %change later:
                sx=loc.locprecnm(sortind);
                if ~isempty(loc.locprecznm)
                sy=loc.locprecznm(sortind);
                else
                    sy=sx;
                end
                loc.sx=sx;
                loc.sy=sy;
                loc.y=yrot(sortind)+zmean;
                loc.znm=loc.znm(sortind);
                for kc=1:length(renderfield)
                    if ~isempty(loc.(renderfield{kc}))
                        loc.(renderfield{kc})=loc.(renderfield{kc})(sortind);
                    end
                end
        end

        end
        
        function rotate_callback(obj,button,b)
            global SMAP_stopnow
            bh=obj.guihandles.rotateb;
            if bh.Value
                bh.FontWeight='bold';
            else
                bh.FontWeight='normal';
            end
             p=obj.getGuiParameters;
            switch p.raxis.selection
                case 'vertical'
                    roih=obj.getPar('sr_roihandle');
                    
                    while bh.Value &&~SMAP_stopnow && strcmp(p.raxis.selection,obj.getSingleGuiParameter('raxis').selection)
                        pos=roih.getPosition;
                        posr=rotpos(pos,obj.getSingleGuiParameter('dangle')*pi/180);
                        roih.setPosition(posr);

                    end           
                case 'horizontal'
                    while bh.Value &&~SMAP_stopnow && strcmp(p.raxis.selection,obj.getSingleGuiParameter('raxis').selection)
                        theta=obj.getSingleGuiParameter('theta');
                        theta=theta-obj.getSingleGuiParameter('dangle')*pi/180;
                        theta=mod(theta,2*pi);                       
                        obj.setGuiParameters(struct('theta',theta));
                        obj.redraw;
                    end   
            end
            
            if  ~ strcmp(p.raxis.selection,obj.getSingleGuiParameter('raxis').selection)
                rotate_callback(obj,button,b)
            end            
        end
        function savemovie_callback(obj,a,b)
            global SMAP_stopnow
            [path,fo]=fileparts(obj.locData.files.file(1).name);
            [file,path]=uiputfile([path filesep fo '.tif']);
            if ~file
                return
            end
            p=obj.getGuiParameters(false,true);
            if length(p.anglerange)==1
                p.anglerange(2)=p.anglerange(1);
                p.anglerange(1)=0;
            end
            if p.dangle<0
                angles=(p.anglerange(2):p.dangle:p.anglerange(1))*pi/180;
            else
                angles=(p.anglerange(1):p.dangle:p.anglerange(2))*pi/180;
            end
            s=size(obj.currentimage.image);
            if length(s)==2
                s(3)=1;
            end
            outim=zeros(s(1),s(2),s(3),length(angles));
            obj.redraw;
            for k=1:length(obj.currentimage.imax)
                obj.recpar{k}.imaxtoggle=false;
                obj.recpar{k}.imax=obj.currentimage.imax(k);
            end

            switch p.raxis.selection
                case 'vertical'
                    roih=obj.getPar('sr_roihandle');
                    pos=roih.getPosition;
                    for k=1:length(angles) 
                        posr=rotpos(pos,angles(k));
                        roih.setPosition(posr);
                        outim(:,:,:,k)=obj.currentimage.image;
                        drawnow
                        if SMAP_stopnow
                            break
                        end
                    end

                case 'horizontal'  
                    for k=1:length(angles) 
                        obj.setGuiParameters(struct('theta',angles(k)));
                        obj.redraw;
                        outim(:,:,:,k)=obj.currentimage.image;
                        drawnow
                        if SMAP_stopnow
                            break
                        end
                    end    
            end
            options.color=true;
            options.message=true;
            options.comp='lzw';

            imout=uint8(outim*(2^8-1));
            saveastiff(imout,[path,file],options)
            obj.recpar={};
        end
        
        function resetazimuth(obj, a,b)
            obj.setGuiParameters(struct('theta',0));
            obj.redraw;
        end
                   
    end
end

function pos=rotpos(pos,angle)
roivec=pos(2,:)-pos(1,:);
mpos=mean(pos,1);
[dx,dy]=rotcoord(roivec(1)/2,roivec(2)/2,angle);
pos(1,1)=mpos(1)-dx;
pos(2,1)=mpos(1)+dx;
pos(1,2)=mpos(2)-dy;
pos(2,2)=mpos(2)+dy;
end

function p=getstereosettings(p,channel)
switch p.stereo.Value
    case 2
        if channel==1
            p.lut.selection='cyan';
        else
            p.lut.selection='red';
        end
        
    case 3
    case 4
end
end

function pard=guidef(obj)
% pard.text1.object=struct('String','parameters','Style','text');
% pard.text1.position=[1,1];

pard.text2.object=struct('String','zmin','Style','text');
pard.text2.position=[1,1];
pard.text3.object=struct('String','zmax','Style','text');
pard.text3.position=[2,1];
pard.setpixelsize.object=struct('String','set pixelsize (x z): ','Style','checkbox','Value',1);
pard.setpixelsize.position=[4,1];
pard.setpixelsize.Width=1.5;
pard.transparencymode.object=struct('String',{{'maximum intensity', 'transparency','balls'}} ,'Style','popupmenu');
pard.transparencymode.position=[6,1];
pard.transparencymode.Width=1.5;
pard.thetat.object=struct('String','Azimuth angle theta ','Style','pushbutton','Callback',@obj.resetazimuth);
pard.thetat.position=[3,1];
pard.thetat.Width=1.1;

pard.zmin.object=struct('Style','edit','String','-400'); 
pard.zmin.position=[1,2.1];
pard.zmin.Width=0.5;
pard.zmax.object=struct('Style','edit','String','400'); 
pard.zmax.position=[2,2.1];
pard.zmax.Width=0.5;
pard.theta.object=struct('Style','edit','String','0'); 
pard.theta.position=[3,2.1];
pard.theta.Width=0.5;
pard.transparencypar.object=struct('Style','edit','String','1'); 
pard.transparencypar.position=[6,2.5];
pard.transparencypar.Width=0.5;
pard.pixrecset.object=struct('Style','edit','String','5 5'); 
pard.pixrecset.position=[4,2.1];
pard.pixrecset.Width=0.5;

pard.showcontrols.object=struct('String','Show Controls','Style','pushbutton','Callback',@obj.showpanel_callback);
pard.showcontrols.position=[1,3];
pard.showcontrols.Width=2;
% pard.ttranslation.object=struct('String','translate','Style','text');
% pard.ttranslation.position=[4,3];
% 
% pard.trot.object=struct('String','rot (command)','Style','text');
% pard.trot.position=[4,4];
% 
% pard.tzoom.object=struct('String','zoom (alt)','Style','text');
% pard.tzoom.position=[8,4];

pard.rotateb.object=struct('String','Rotate','Style','togglebutton','Callback',@obj.rotate_callback);
pard.rotateb.position=[2,3];
pard.raxis.object=struct('String',{{'horizontal','vertical'}},'Style','popupmenu');%,'Callback',@obj.axischange_callback);
pard.raxis.position=[3,3];
pard.danglet.object=struct('String','step','Style','text');
pard.danglet.position=[3,4];
pard.danglet.Width=0.5;
pard.dangle.object=struct('String','3','Style','edit');
pard.dangle.position=[3,4.5];
pard.dangle.Width=0.5;

pard.savemovie.object=struct('String','save movie','Style','pushbutton','Callback',@obj.savemovie_callback);
pard.savemovie.position=[4,3];
pard.tx.object=struct('String','min max angle (deg)','Style','text');
pard.tx.position=[5,3];
pard.anglerange.object=struct('String','0 360','Style','edit');
pard.anglerange.position=[5,4];

pard.stereo.object=struct('Style','popupmenu','String',{{'no stereo','anaglyphs (color)','free-view','cross-eyed'}}); 
pard.stereo.position=[7,1];
pard.stereo.Width=1;
pard.tx2.object=struct('String','pixel (mm)','Style','text');
pard.tx2.position=[7,2];
pard.tx2.Width=.6;
pard.monitorpmm.object=struct('Style','edit','String','0.3'); 
pard.monitorpmm.position=[7,2.6];
pard.monitorpmm.Width=.4;

pard.tx3.object=struct('String','d plane','Style','text');
pard.tx3.position=[7,3];
pard.tx3.Width=.6;
pard.dplanemm.object=struct('Style','edit','String','1000');  
pard.dplanemm.position=[7,3.6];
pard.dplanemm.Width=.4;

pard.tx4.object=struct('String','M two images','Style','text');
pard.tx4.position=[7,4];
pard.tx4.Width=.6;
pard.stereomagnification.object=struct('Style','edit','String','1');  
pard.stereomagnification.position=[7,4.6];
pard.stereomagnification.Width=.4;
end