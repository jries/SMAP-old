classdef VersatileRenderer<interfaces.DialogProcessor
    methods
        function obj=VersatileRenderer(varargin)        
            obj@interfaces.DialogProcessor(varargin{:}) ;
            obj.inputParameters={'sr_layerson','numberOfLayers'};
        end
        
        function out=run(obj,p)
            out=[];
            ax=initaxis(p.resultstabgroup,'image');
            lochere=obj.locData.copy;
            [~,inroi]=lochere.getloc('xnm','position','roi');
            lochere.removelocs(~inroi);
            lochere.regroup;
            v1=lochere.loc.(p.assignfield1.selection)/p.pixelsize1;
            v2=lochere.loc.(p.assignfield2.selection)/p.pixelsize2;
            lochere.loc.xnm=v1;
            lochere.loc.ynm=v2;

            v1=lochere.grouploc.(p.assignfield1.selection)/p.pixelsize1;
            v2=lochere.grouploc.(p.assignfield2.selection)/p.pixelsize2;
            lochere.grouploc.xnm=v1;
            lochere.grouploc.ynm=v2;
            
            pall=obj.getLayerParameters;
            phere.gaussfac=0;
            phere.mingausspix=p.filtersize;
            phere.sr_pixrec=1;

            phere.rangex=[p.min1 p.max1]/p.pixelsize1;
            phere.rangey=[p.min2 p.max2]/p.pixelsize2;
            phere.mingausspix=p.filtersize;
            phere.gaussfac=.01;
            phere.mingaussnm=0;
            
            phere.sr_axes=ax;
            for k=1:length(pall)
                pall{k}=copyfields(pall{k},p);
                pall{k}=copyfields(pall{k},phere);
            end
            oovim=TotalRender(lochere,pall,{'xnm','ynm'});
            
%             srim=renderAnyField(obj.locData,p);
%             
        end
        function pard=guidef(obj)
            pard=guidef;
        end
        function initGui(obj)
            obj.guihandles.assignfield1.Callback={@assignfield_callback,obj,1};
             obj.guihandles.assignfield2.Callback={@assignfield_callback,obj,2};
            callobj=obj;
            obj.addSynchronization('locFields',obj.guihandles.assignfield1,'String');%,{@callobj.assignfield_callback,1});
            obj.addSynchronization('locFields',obj.guihandles.assignfield2,'String');%,{@callobj.assignfield_callback,2});
        end
    end
end
function assignfield_callback(object,a,obj,axis)

%     object=obj.guihandles.(['assignfield' int2str(axis)]);
    field=object.String{object.Value};
    v=obj.locData.getloc(field,'layer',1,'position','roi').(field);
    q=myquantile(v,[0.0001,0.5,0.9999]);

    dx=(q(3)-q(1))/100;
    dx=10^floor(log10((dx)));
    obj.guihandles.(['pixelsize' num2str(axis)]).String=num2str(dx);
    minv=floor(q(1)/dx)*dx-5*abs(dx);
    maxv=ceil(q(3)/dx)*dx+5*abs(dx);
    obj.guihandles.(['min' num2str(axis)]).String=num2str(minv);
    obj.guihandles.(['max' num2str(axis)]).String=num2str(maxv);

% changemode(0,0,obj);
%setSingeleParR
end
function pard=guidef
pard.t1.object=struct('String','horizontal axis','Style','text');
pard.t1.position=[2,2];
pard.t2.object=struct('String','vertical axis','Style','text');
pard.t2.position=[2,3];

pard.assignfield1.object=struct('Style','popupmenu','String',{{'n1,n2'}});
pard.assignfield1.position=[3,2];
pard.assignfield2.object=struct('Style','popupmenu','String',{{'n1,n2'}});
pard.assignfield2.position=[3,3];
pard.assignfield1.object.TooltipString='choose which field to render';
pard.assignfield2.object.TooltipString=pard.assignfield1.object.TooltipString;

pard.t3.object=struct('String','pixelsize','Style','text');
pard.t3.position=[4,1];
pard.t4.object=struct('String','min','Style','text');
pard.t4.position=[5,1];
pard.t5.object=struct('String','max','Style','text');
pard.t5.position=[6,1];

pard.pixelsize1.object=struct('String','','Style','edit');
pard.pixelsize1.position=[4,2];
pard.pixelsize2.object=struct('String','','Style','edit');
pard.pixelsize2.position=[4,3];

pard.min1.object=struct('String','','Style','edit');
pard.min1.position=[5,2];
pard.min2.object=struct('String','','Style','edit');
pard.min2.position=[5,3];

pard.max1.object=struct('String','','Style','edit');
pard.max1.position=[6,2];
pard.max2.object=struct('String','','Style','edit');
pard.max2.position=[6,3];

pard.t6.object=struct('String','sigma rec (pix)','Style','text');
pard.t6.position=[2,4];
pard.filtersize.object=struct('String','1','Style','edit');
pard.filtersize.position=[3,4];
end

