function [locsout,indcombined,hroio]=getlocs(locData,fields,varargin)
%returns subset of localizations. 
%[locs,handle to roi]=getloc(obj,{fields},'PropertyName','Property')
%fields: any locData field, in addition: 'ingrouped',
%'inungrouped','within'
% inlayeru, inlayerg: cell array with indices for each layer

%Properties:
%'grouping': 'grouped', 'ungrouped' (second default),'layer'(default). If not set and combined with
%'layer, use individual layer grouping. If any of layers is ungrouped,
%return ungrouped localizations.
%'channels': double vector of channels
%'filenumber': double vector of filenumbers
%'position': 'all' (default),'roi','fov': position filter
%'position': vector: [centerx, centery, widhtx widthy]
%'layer': double number or vector of layers.

%'removeFilter': cell array of filter names to remove 
%'within', indices which localizations to consider.

 p=roiparser(varargin); 
if ischar(fields)
    p.fields={fields};
else
    p.fields=fields;
end
 indlayer={};
if isempty(locData.loc)
    for k=1:length(p.fields)
        locsout.(p.fields{k})=[];
    end
    indcombined=[];hroio=[];
    return
end


 %check if empty. necessary at beginning where xnm etc can be empty.
if isempty(p.channel)&&strcmpi(p.grouping,'layer')&&isempty(p.layer)&&strcmpi(p.position,'default')&&isempty(p.filenumber)&&isempty(p.removeFilter)&&isempty(p.within)
locs=locData.loc;
    for k=1:length(p.fields)
        if isfield(locs,p.fields{k})
            locsout.(p.fields{k})=locs.(p.fields{k});
        elseif strcmp(p.fields{k},'ingrouped')
            locsout.(p.fields{k})=true(length(locData.grouploc.frame),1);
        elseif strcmp(p.fields{k},'inungrouped')
            locsout.(p.fields{k})=true(length(locData.loc.frame),1);
        else
            locsout.(p.fields{k})=[];
        end
    end
    hroio=[];
    indcombined=[];
    return
end

 hroio=[];
%  fields=fieldnames(locData.loc);           
         
%  if ~isempty(p.layer)  %determine group from layer   
%       %determine output format: only if all grouped, return grouped data.
%      if strcmp(p.grouping, 'layer')
%          groupchecked=true;       
%          for k=1:length(p.layer)
%             group=group&locData.getPar(['layer' num2str(p.layer(k)) '_groupcheck']);
%          end
%          if group
%              p.grouping='grouped';
%              grouping=true;
%          else
%              p.grouping='ungrouped';
%              grouping=false;
%          end
%      end
%  end
 
 switch p.grouping
     case {'ungrouped',false}
         grouping=false;
         locs=locData.loc;
     case {'grouped',true}
         if isempty(locData.grouploc)
             locData.regroup;
         end
         locs=locData.grouploc;
         grouping=true;
     case 'layer'
          if isempty(p.layer)
              grouping=false;
              locs=locData.loc;
          else
               group=true;       
                 for k=1:length(p.layer)
                     grouplayer(k)=locData.getPar(['layer' num2str(p.layer(k)) '_groupcheck']);
                    group=group&grouplayer(k);
                 end
                 if group
                     p.grouping='grouped';
                     grouping=true;
                     locs=locData.grouploc;
                 else
                     p.grouping='ungrouped';
                     grouping=false;
                     locs=locData.loc;
                 end
          end
         
 end

 %get filters for layers
% remove filter in case field is defined

indfilter=false(size(locs.xnm));
% layerused=false;
for k=1:length(p.layer)
    if p.layer(k)>length(locData.layer)||p.layer(k)<1
        disp('getloc layer out of range')
        continue
    end
     filterold=locData.getFilter(p.layer(k),grouping); %use localizations if in any of the layers
     if ~iscell(p.removeFilter)
         p.removeFilter={p.removeFilter};
     end
      if isnumeric(p.position)||strcmp(p.position,'all')||strcmp(p.position,'default')
          p.removeFilter=[p.removeFilter,{'xnm','ynm'}];
      end
      if ~isempty(p.channel)
          p.removeFilter=[p.removeFilter,{'channel'}];
      end
      if ~isempty(p.filenumber)
          p.removeFilter=[p.removeFilter,{'filenumber'}];
      end
      if ~isempty(p.removeFilter)
         filter=myrmfield(filterold,p.removeFilter);
         locData.setFilter(filter,p.layer(k),grouping);
      end
      indlayer{k}=locData.inFilter(p.layer(k),grouping);
%      indfilter=indfilter|locData.inFilter(p.layer(k),grouping);
     indfilter=indfilter|indlayer{k};
     locData.setFilter(filterold,p.layer(k),grouping);
%      layerused=true;
end 
 if isempty(p.layer)
     indfilter=true(size(locs.xnm));
 end

 if ~isempty(p.channel)
     indfilterh=(locs.channel==p.channel(1));
     for k=2:length(p.channel)
        indfilterh=indfilterh|(locs.channel==p.channel(k));
     end
     indfilter=indfilter&indfilterh;
 end
 if ~isempty(p.filenumber)
     indfilterh=(locs.filenumber==p.filenumber(1));
     for k=2:length(p.filenumber)
        indfilterh=indfilterh|(locs.filenumber==p.filenumber(k));
     end
     indfilter=indfilter&indfilterh;
 end
     
if isnumeric(p.position)
     pos=p.position(1:2);
        sr_size=p.position(3:4)/2;
        indpos=locs.xnm>pos(1)-sr_size(1) & locs.xnm<pos(1)+sr_size(1) & locs.ynm>pos(2)-sr_size(2) & locs.ynm<pos(2)+sr_size(2);
else
    switch p.position
        case {'all','default'}
            indpos=true(size(locs.xnm));
        case 'roi'
            [indpos,hroio,strucout]=getinroi(locData,locs.xnm,locs.ynm);
            if isfield(strucout,'xnmline')
%                 locs.xnmline=zeros(size(locs.xnm));
%                 locs.ynmline=zeros(size(locs.ynm));
%                 locs.xnmline(indpos)=strucout.xnmline;
%                 locs.ynmline(indpos)=strucout.ynmline;
               locs.xnmline=strucout.xnmline;
               locs.ynmline=strucout.ynmline;
            end
        case 'fov'
            pos=locData.getPar('sr_pos');
            sr_size=locData.getPar('sr_size');
            indpos=locs.xnm>pos(1)-sr_size(1) & locs.xnm<pos(1)+sr_size(1) & locs.ynm>pos(2)-sr_size(2) & locs.ynm<pos(2)+sr_size(2);
        otherwise %numerical position vector
           disp('getlocs: position description not valid')
    end
end

if ~isempty(p.within)
    indwithin=p.within;
    if length(indwithin)<length(indfilter)
        indwithin=grouped2ungrouped(locData,indwithin);
    elseif length(indwithin)>length(indfilter)
        indwithin=ungrouped2grouped(locData,indwithin);
    end
else
    indwithin=1;
end

indcombined=indfilter&indpos&indwithin;

if isempty(p.fields)||any(strcmpi(p.fields,'all'))
    p.fields=fieldnames(locs);
end

 for k=1:length(p.fields)
     field=p.fields{k};
     if isfield(locs,field)
        locsout.(field)=locs.(field)(indcombined);
     elseif strcmp(p.fields{k},'ingrouped')
            locsout.(p.fields{k})=getindices(locData,indcombined,1);
     elseif strcmp(p.fields{k},'inungrouped')
            locsout.(p.fields{k})=getindices(locData,indcombined,0);
     elseif strcmp(p.fields{k},'inlayeru')
            for l=1:length(indlayer)
                indlayer{l}=getindices(locData,indcombined,0)&getindices(locData,indlayer{l} ,0);
            end
            locsout.(p.fields{k})=indlayer;
     elseif strcmp(p.fields{k},'inlayerg')
            for l=1:length(indlayer)
                indlayer{l}=getindices(locData,indcombined,1)&getindices(locData,indlayer{l},1);
            end
            locsout.(p.fields{k})=indlayer;            
     else
         locsout.(field)=[];
     end
 end            

 
end
function ind=getindices(obj,indcombined,isgrouped)
    if isgrouped
        fn=fieldnames(obj.grouploc);
        len=length(obj.grouploc.(fn{1}));
    else
        fn=fieldnames(obj.loc);
        len=length(obj.loc.(fn{1}));
    end
    if len==length(indcombined)
        ind=indcombined;
    elseif isgrouped
        ind=ungrouped2grouped(obj,indcombined);
    else
        ind=grouped2ungrouped(obj,indcombined);
    end   
end

function ind=ungrouped2grouped(obj,indcombined)
gind=unique(obj.loc.groupindex(indcombined));
ind=false(length(obj.grouploc.frame),1);
ind(gind)=true;
end

function ind=grouped2ungrouped(obj,indcombined)
ind=indcombined(obj.loc.groupindex);
end

function [indg,hroio,strucout]=getinroi(obj,vx,vy)
%obj=locData
hroi=obj.getPar('sr_roihandle');
% hroi=obj.parameters.roihandle;
strucout=[];
if isa(hroi,'imroi')&&isvalid(hroi)  
        
    if isa(hroi,'imline')
        lw=obj.getPar('linewidth_roi')/2;
        pol=hroi.getPosition;
        dpol=pol(2,:)-pol(1,:);
        alpha=-atan2(dpol(1),dpol(2));

        len=sqrt(sum(dpol.^2))*1000/2;
        midp=mean(pol,1)*1000;
        [xr,yr]=rotcoord(vx-midp(1),vy-midp(2),alpha);
%         xs=vx-midp(1);ys=vy-midp(2);
%         xr=xs*cos(alpha)+ys*sin(alpha);
%         yr=ys*cos(alpha)-xs*sin(alpha);
%         indg=abs(xr)<lw&abs(yr)<len;
        indb=abs(xr)>lw|abs(yr)>len;
        indg=~indb;
        xr(indb)=0;yr(indb)=0;
%         xr=xr.*indg;
%         yr=yr.*indg;
        
%         strucout.xnmline=zeros(length(indg),1);
%         strucout.ynmline=zeros(length(indg),1);
        strucout.ynmline=xr;strucout.xnmline=yr;
%  strucout.ynmline(indg)=xr(indg);strucout.xnmline(indg)=yr(indg);
    elseif isa(hroi,'impoint')
        lw=obj.parameters.linewidth_roi/2;
        pol=hroi.getPosition*1000;
        indg=abs(vx-pol(1))<lw&abs(vy-pol(2))<lw;
    else
        imbw=createMask(hroi,obj.getPar('sr_imagehandle'));
        sizeim=size(imbw);
        pos=obj.getPar('sr_pos');
        sizesr=obj.getPar('sr_size');
        pixrec=obj.getPar('sr_pixrec');
        xl=min(max(1,round((vx-pos(1)+sizesr(1))/pixrec)),sizeim(2));
        yl=min(max(1,round((vy-pos(2)+sizesr(2))/pixrec)),sizeim(1));
        linind=sub2ind(sizeim,yl,xl);
        indg=imbw(linind);
    end
    hroio=hroi;
else 
    srp=obj.getPar('sr_pos');
    ssize=obj.getPar('sr_size');
    pos=(srp(1:2)-ssize);
    
    hroio.getPosition=[pos(1) pos(2) 2*ssize(1) 2*ssize(2)]/1000;
    indg=vx>pos(1)&vx<pos(1)+2*ssize(1)&vy>pos(2)&vy<pos(2)+2*ssize(2);
    strucout=[];
end
end

function pres=roiparser(args)
% fields{end+1}='all';
p = inputParser;   
p.KeepUnmatched=true;
% addOptional(p,'fields',{},@(x)  any(myvalidatestring(x,fields)));
addParameter(p,'grouping','layer',@(x) any(myvalidatestring(x,{'grouped','ungrouped','layer'})));
addParameter(p,'layer',[],@isnumeric);
addParameter(p,'channel',[],@isnumeric);
addParameter(p,'filenumber',[],@isnumeric);
addParameter(p,'position','default');
addParameter(p,'removeFilter',{});
addParameter(p,'within',[]);
parse(p,args{:});
pres=p.Results;
end