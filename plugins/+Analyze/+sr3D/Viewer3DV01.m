classdef Viewer3DV01<interfaces.DialogProcessor
    properties
        axis
        timer
        theta=0;
%         locCopy;
    end
    methods
        function obj=Viewer3DV01(varargin)        
            obj@interfaces.DialogProcessor(varargin{:}) ;
            obj.inputParameters=[renderSMAP drawerSMAP displayerSMAP];
            obj.inputParameters{end+1}='sr_roihandle';
            obj.inputParameters{end+1}='linewidth_roi';
            obj.inputParameters{end+1}='layers';
            obj.inputParameters{end+1}='numberOfLayers';
            
             obj.showresults=true;

        end
        function makeGui(obj)
            makeGui@interfaces.DialogProcessor(obj);
            h=obj.guihandles;
            
            h.ptranslation=makepanel(h.ttranslation,'translation','');
            h.protation=makepanel(h.trot,'rot (command/strg)','command');
            h.pzoom=makepanel(h.tzoom,'zoom (alt)','alt');
%             h.protation=uipanel('Parent',obj.handle,'Units','pixels','Position',pos,'Title','translation');
%             h.pzoom=uipanel('Parent',obj.handle,'Units','pixels','Position',pos,'Title','translation');
            
            function h=makepanel(htext,title,modifier)
                pos=htext.Position;
                pos(4)=pos(4)*4;
                h=uipanel('Parent',obj.handle,'Units','pixels','Position',pos,'Title',title);
                h.Units='normalized';
                uicontrol('Parent',h,'Units','normalized','Position',[1 1 1 1]/3,'String','0','Callback',{@obj.keypress,struct('Modifier',modifier,'Key','0')})
                uicontrol('Parent',h,'Units','normalized','Position',[0 1 1 1]/3,'String','<-','Callback',{@obj.keypress,struct('Modifier',modifier,'Key','leftarrow')})
                uicontrol('Parent',h,'Units','normalized','Position',[2 1 1 1]/3,'String','->','Callback',{@obj.keypress,struct('Modifier',modifier,'Key','rightarrow')})
                uicontrol('Parent',h,'Units','normalized','Position',[1 2 1 1]/3,'String','^','Callback',{@obj.keypress,struct('Modifier',modifier,'Key','uparrow')})
                uicontrol('Parent',h,'Units','normalized','Position',[1 0 1 1]/3,'String','v','Callback',{@obj.keypress,struct('Modifier',modifier,'Key','downarrow')})   
                uicontrol('Parent',h,'Units','normalized','Position',[0 2 1 1]/3,'String','x (,)','Callback',{@obj.keypress,struct('Modifier',modifier,'Key','comma')})
                uicontrol('Parent',h,'Units','normalized','Position',[2 2 1 1]/3,'String','(.)','Callback',{@obj.keypress,struct('Modifier',modifier,'Key','period')})                  
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
%             obj.axis=initaxis(obj.resultstabgroup,'3Dimage');
            obj.axis.Units='normalized';
            obj.axis.Position=[0.05 0.05 .95 .9];
            
             set(obj.axis,'NextPlot','replacechildren','PickableParts','all','Units','pixels')
            fig=getParentFigure(obj.axis);
             axis(obj.axis,'tight')
             axis(obj.axis,'equal')
             axis(obj.axis,'ij')
             
             set(fig,'WindowKeyPressFcn',{@obj.keypress,[]})
             obj.theta=0;
             obj.timer=uint64(0);
             obj.redraw
             

        end

        function pard=pardef(obj)
            pard=pardef;
        end
        function keypress(obj,a,d2,data)
            if isempty(data)
                data=d2;
            end
            
%             tic
            %no modifier:move
            %commamd: rotate
            %control: thickness
            %shift: small step
%             data.Key
%  data.Key
            p=obj.getGuiParameters;
            if (length(data.Key)<5 && ~strcmp(data.Key,'0'))||(length(data.Key)>5 && ~(strcmp(data.Key(end-4:end),'arrow')|| strcmp(data.Key,'period')|| strcmp(data.Key,'comma')))
                return
            end
           
            roih=obj.getPar('sr_roihandle');
            pos=roih.getPosition;
            roivec=pos(2,:)-pos(1,:);
            roivecp(2)=roivec(1);
            roivecp(1)=-roivec(2);
            step=0.1;
            stepl=0.3;
            dphi=pi/32;
            dtheta=pi/8;
            if any(strcmp(data.Modifier,'shift'))
                stepfac=0.2;
            else
                stepfac=1;
            end
%             data.Modifier
            if any(strcmp(data.Modifier,'command'))||any(strcmp(data.Modifier,'control'))
                %rotate
                phi=0;
                switch data.Key
                    
                    case 'uparrow'
                       %tilt up down
                       obj.theta=obj.theta+dtheta*stepfac;
                       if obj.theta>pi
                           obj.theta=obj.theta-2*pi;
                       end
                       if obj.theta<-pi
                           obj.theta=obj.theta+2*pi;
                       end
                       phi=0;
                    case 'downarrow'
                        phi=0;
                       obj.theta=obj.theta-dtheta*stepfac;
                       if obj.theta>pi
                           obj.theta=obj.theta-2*pi;
                       end
                       
                       if obj.theta<-pi
                           obj.theta=obj.theta+2*pi;
                       end
                    case 'leftarrow'
                       phi=dphi*stepfac;
                    case 'rightarrow'
                        phi=-dphi*stepfac;
                    case '0'
                        obj.theta=0;
                end
                 mpos=mean(pos,1);
                [dx,dy]=rotcoord(roivec(1)/2,roivec(2)/2,phi);
                pos(1,1)=mpos(1)-dx;
                pos(2,1)=mpos(1)+dx;
                pos(1,2)=mpos(2)-dy;
                pos(2,2)=mpos(2)+dy;
            elseif any(strcmp(data.Modifier,'alt'))
                %change size
                switch data.Key
                    case 'period'
                        lw=obj.getPar('linewidth_roi');
                        lw2=lw*(1+step*stepfac);
                        obj.setPar('linewidth_roi',lw2);
                       %tilt up down
                    case 'comma'
                        lw=obj.getPar('linewidth_roi');
                        lw2=lw*(1-step*stepfac);
                        obj.setPar('linewidth_roi',lw2);
                      
                    case 'leftarrow'
                        
                       pos(1,:)=pos(1,:)+roivec/2*step*stepfac;
                       pos(2,:)=pos(2,:)-roivec/2*step*stepfac;
                    case 'rightarrow'
                        pos(1,:)=pos(1,:)-roivec/2*step*stepfac;
                       pos(2,:)=pos(2,:)+roivec/2*step*stepfac;
                    case 'uparrow'
                        po.zmin=p.zmin-stepfac*step*(p.zmax-p.zmin);
                        po.zmax=p.zmax+stepfac*step*(p.zmax-p.zmin);
                        obj.setGuiParameters(po);
                    case 'downarrow'
                        po.zmin=p.zmin+stepfac*step*(p.zmax-p.zmin);
                        po.zmax=p.zmax-stepfac*step*(p.zmax-p.zmin);
                        obj.setGuiParameters(po);
                    
                       
                end
            else
                switch data.Key
                    case 'period'
                        lw=obj.getPar('linewidth_roi')/1000;
                        pos(1,:)=pos(1,:)+stepl*roivecp./norm(roivecp)*stepfac*lw;
                        pos(2,:)=pos(2,:)+stepl*roivecp./norm(roivecp)*stepfac*lw;
                    case 'comma'
                        lw=obj.getPar('linewidth_roi')/1000;
                        pos(1,:)=pos(1,:)-stepl*roivecp./norm(roivecp)*stepfac*lw;
                        pos(2,:)=pos(2,:)-stepl*roivecp./norm(roivecp)*stepfac*lw;
                    case 'leftarrow'
                        pos(1,:)=pos(1,:)+step*roivec*stepfac;
                        pos(2,:)=pos(2,:)+step*roivec*stepfac;
                    case 'rightarrow'
                        pos(1,:)=pos(1,:)-step*roivec*stepfac;
                        pos(2,:)=pos(2,:)-step*roivec*stepfac;
                    case 'uparrow'
                        po.zmin=p.zmin+stepfac*step*(p.zmax-p.zmin);
                        po.zmax=p.zmax+stepfac*step*(p.zmax-p.zmin);
                        obj.setGuiParameters(po);
                    case 'downarrow'
                        po.zmin=p.zmin-stepfac*step*(p.zmax-p.zmin);
                        po.zmax=p.zmax-stepfac*step*(p.zmax-p.zmin);
                        obj.setGuiParameters(po);
                    case '0'
                        po.zmin=p.zmin-(p.zmax+p.zmin)/2;
                        po.zmax=p.zmax-(p.zmax+p.zmin)/2;
                        obj.setGuiParameters(po);
                end
            end
            
            
            roih.setPosition(pos);
%             toc
        end
        function redraw(obj)
            
            if toc(obj.timer)<0.01
                return
            end
             p=obj.getAllParameters;
%             tic
            locCopy=obj.locData; %maybe not needed
            lo=logical(obj.getPar('sr_layerson'));
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
            if group(1)
                [loc,indu]=locCopy.getloc({'xnmline','ynmline','znm','locprecnm','locprecznm',renderfield{:}},'position','roi','grouping','ungrouped');
                
                [yrot,depth]=rotcoord(loc.znm-zmean,loc.ynmline,obj.theta);
%                 [zmrot]=rotcoord(zmean,0,obj.theta);
                depth=depth-min(depth);
%                 md=max(depth);
                [~,sortind]=sort(depth);
                loc.x=loc.xnmline(sortind);
                %change later:
                sx=loc.locprecnm(sortind);
                sy=loc.locprecznm(sortind);
                loc.sx=sx(sortind);
                loc.sy=sy(sortind);
                loc.y=yrot(sortind)+zmean;
                loc.znm=loc.znm(sortind);
                for k=1:length(renderfield)
                    if ~isempty(loc.(renderfield{k}))
                        loc.(renderfield{k})=loc.(renderfield{k})(sortind);
                    end
                end
    %             locCopy.loc.x=locCopy.loc.xnmline;
%                 loc.intensity_render=intd(sortind);
                
            end
            
            if group(2)
                [locg,indg]=locCopy.getloc({'xnmline','ynmline','znm','locprecnm','locprecznm',renderfield{:}},'position','roi','grouping','grouped');
                [yrot,depth]=rotcoord(locg.znm-zmean,locg.ynmline,obj.theta);
    %             yline=zeros(length(indin),1);
    %             yline(indin)=yrot;
                depth=depth-min(depth);
                md=max(depth);
    %             intd=zeros(length(indin),1);
                intd=1./(1+4*depth/md);
                
                [~,sortind]=sort(depth);

                locg.x=locg.xnmline(sortind);
                %change later:
                sx=locg.locprecnm(sortind);
                sy=locg.locprecznm(sortind);
                locg.sx=sx(sortind);
                locg.sy=sy(sortind);
                locg.znm=locg.znm(sortind);
                locg.y=yrot(sortind)+zmean;
                for k=1:length(renderfield)
                    if ~isempty(locg.(renderfield{k}))
                        locg.(renderfield{k})=locg.(renderfield{k})(sortind);
                    end
                end
            end
%             [rz]=rotcoord([p.zmin p.zmax],[-p.linewidth_roi p.linewidth_roi]/2,obj.theta);
            
            ph.rangey=[p.zmin p.zmax];
%             ph.rangey(1)=min(rz);
%             ph.rangey(2)=max(rz);
            
            roih=obj.getPar('sr_roihandle');
            rpos=roih.getPosition;
    %     mpos=mean(rpos,1);
            lr=sqrt(sum((rpos(2,:)-rpos(1,:)).^2));
            rx=[-lr/2 lr/2]*1000;
%             ry=[-0.5 0.5]*p.linewidth_roi;
        
            ax=obj.axis;
%             disp('redraw callback')
            obj.timer=tic;
           
            ph.sr_roihandle=obj.getPar('sr_roihandle');
%             p.rangey=[-300 300];
            
            ph.rangex=rx;
            if p.settransparency
            transparency=p.transparency;
            else
                transparency=[];
            end
            ph.sr_axes=[];
            if p.setpixelsize
                ph.sr_pixrec=p.pixrecset;
            end
            for k=1:p.numberOfLayers
                pl=p.(['layer' num2str(k) '_']);

                if pl.layercheck

%                      pl.mingaussnm=0;%why????
                     pr=copyfields(copyfields(p,pl),ph);
                     if pl.groupcheck
                         indroi=obj.locData.getloc('ingrouped','layer',k,'position','roi').ingrouped;  
                        layer(k).images.srimage=renderSMAP(locg,pr,k,indroi(indg),transparency);
                     else
                         indroi=obj.locData.getloc('inungrouped','layer',k,'position','roi').inungrouped;  
                         layer(k).images.srimage=renderSMAP(loc,pr,k,indroi(indu),transparency);
                     end

%                     layer(k).images.srimage.rangex=layer(k).images.srimage.rangex;
%                     layer(k).images.srimage.rangey=layer(k).images.srimage.rangey;
                    layer(k).images.finalImages=drawerSMAP(layer(k).images.srimage,pr);        

                end
            end
            srim=displayerSMAP(layer,pr);
            
%             [zim]=anyRender(locCopy,p,'x','xnmline','y','ynmrot','sx','locprecnm','sy','locprecznm','within',indin,'position','roi','groupstate',group);
            imagesc(ph.rangex,ph.rangey,srim.image,'Parent',ax);
            title(obj.theta,'Parent',ax)
            drawnow
%             toc

    %directly call renderSMAP, circumvent any_render
%             toc(obj.timer)
        end
    end
end


function pard=pardef
pard.text1.object=struct('String','parameters','Style','text');
pard.text1.position=[1,1];

pard.text2.object=struct('String','zmin','Style','text');
pard.text2.position=[2,1];
pard.text3.object=struct('String','zmax','Style','text');
pard.text3.position=[3,1];
pard.setpixelsize.object=struct('String','set pixelsize (x z): ','Style','checkbox','Value',1);
pard.setpixelsize.position=[4,1];
pard.setpixelsize.Width=1.5;
pard.settransparency.object=struct('String','use transparency: ','Style','checkbox');
pard.settransparency.position=[5,1];
pard.settransparency.Width=1.5;


pard.zmin.object=struct('Style','edit','String','-400'); 
pard.zmin.position=[2,2.1];
pard.zmin.Width=0.5;
pard.zmax.object=struct('Style','edit','String','400'); 
pard.zmax.position=[3,2.1];
pard.zmax.Width=0.5;
pard.transparency.object=struct('Style','edit','String','1'); 
pard.transparency.position=[5,2.1];
pard.transparency.Width=0.5;
pard.pixrecset.object=struct('Style','edit','String','5 5'); 
pard.pixrecset.position=[4,2.1];
pard.pixrecset.Width=0.5;


pard.ttranslation.object=struct('String','translate','Style','text');
pard.ttranslation.position=[4,3];

pard.trot.object=struct('String','rot (command)','Style','text');
pard.trot.position=[4,4];

pard.tzoom.object=struct('String','zoom (alt)','Style','text');
pard.tzoom.position=[8,4];
end