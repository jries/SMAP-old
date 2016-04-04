classdef LocalizationData<interfaces.GuiParameterInterface
    % Stores all localization related information and provides access to
    % localization related functionality
    properties
        loc  %ungrouped localizations. locs.(field) is vector of lenth(number of localizations)
        grouploc %grouped localizations. Not calculated on the fly for performance
        layer %stores layer-specific information = filters.
        files %localization files and file meta data
        SE %siteexplorer object, linked here why?
%         iscopy=true;
    end
    
    methods

        function obj=LocalizationData
            %constructor. calls .clear to initialize 
            obj.clear;
        end
        function clear(obj,part)
            %initializes properties to empty fields, sets files property to
            %no files
            if nargin<2||strcmpi(part,'all')
            obj.loc=[];
            obj.grouploc=[];
             obj.layer=[];
            obj.files.filenumberEnd=0;
            obj.files.file=[];    
            if ~isempty(obj.SE)
                obj.SE.clear;
            end
            elseif strcmpi(part,'filter')
                obj.layer=[];
            end
            obj.inputParameters={'numberOfLayers','filtertable','layer','group_dx','group_dt'};
        end
        function setloc(obj,name, value)
            %Add a field to all localizations. setloc(field, value) sets obj.loc.(name)=value
            obj.loc.(name)=real(value);
        end
        function addloc(obj,name,value)
            %addloc(field, value) Appends localizations to exisiting field, sets localizations
            %if not exisitng. 
            %FIX: can result in field vecotrs of different length. Pad with
            %zeros?
            if isfield(obj.loc,name)
                l1=length(obj.loc.(name));
                l2=length(value);
%                 v1=obj.loc;
%                 value=[obj.loc.(name);real(value)];
                obj.loc.(name)(l1+1:l1+l2)=value;
            else
                obj.setloc(name,value)
            end
        end
        function [locsout,indout,hroi]=getloc(obj,fields,varargin)
            %returns subset of localizations. 
            %[locs,indices of locs,handle to roi]=getloc(obj,{fields},'PropertyName','Property')
            %Properties:
            %'grouping': 'grouped', 'ungrouped' (default)
            %'channels': double vector of channels
            %'filenumber': double vector of filenumbers
            %'position': 'all' (default),'roi','fov': position filter
            %'layer': double number of layer. overrides channel and grouping. Onyl one layer possible
            %(due to grouping)
            %'removeFilter': cell array of filter names to remove 
            [locsout,indout,hroi]=getlocs(obj,fields,varargin{:});
        end
            
        function mask=getRoiMask(obj)
            %returns binary imge of roi: mask=obj.getRoiMask
            hroi=obj.parameters.roihandle;
            mask=createMask(hroi);
        end
        
        function removelocs(obj,ind,removepart)
            %removes localizations with specified indices from all fields
            %obj.removelocs(indices,removepart)
            %optional removepart='loc' (default) or 'grouploc', determine which part to remove
            if nargin<3
                removepart='loc';
            end
            fn=fieldnames(obj.(removepart));
            for k=1:length(fn)
                obj.(removepart).(fn{k})(ind)=[];
            end
        end
        
        function filter=getFilter(obj,layer, grouping)
            %get filter structure from layer, 
            %optional: specify grouping=true or false, if not specified, take from layer
            %filter=obj.getFilter(layer, grouping)
            if nargin<3
                grouping=obj.isgrouped(layer);
            end
            
            if layer<=length(obj.layer)
            if grouping
                filter=obj.layer(layer).groupfilter;
            else
                filter=obj.layer(layer).filter;
            end  
            else
                filter=[];
            end
        end
        
        
        function setFilter(obj,filter,layer,grouping) 
            %set filter structure for specific layer (from getFilter)
            %obj.setFilter(filter,layer,grouping)
            if nargin<4
                grouping=obj.isgrouped(layer);
            end
            if grouping
                obj.layer(layer).groupfilter=filter;
            else
                obj.layer(layer).filter=filter;
            end  
        end
        
        function removeFilter(obj,field,layer)
            %remove specified filter from specified layer
            % removeFilter(field,layer)
            if layer<=length(obj.layer)
            obj.layer(layer).filter=myrmfield(obj.layer(layer).filter,field);
            obj.layer(layer).groupfilter=myrmfield(obj.layer(layer).groupfilter,field);
            end
        end
        
        function indin=inFilter(obj,layer,grouping)
            %returns indices of filtered localizations in specific layer.
            %Grouping (t/f) can be specified or is taken from layer
            %indices=inFilter(layer,grouping)
            if nargin<3
                grouping=obj.isgrouped(layer);
            end
            
            if grouping
                fn=fieldnames(obj.grouploc);
                indin=true(length(obj.grouploc.(fn{1})),1);
            else
                fn=fieldnames(obj.loc);
                indin=true(length(obj.loc.(fn{1})),1);
            end
            filters=obj.getFilter(layer,grouping); 
            if ~isempty(filters)
            fields=fieldnames(filters);
                if ~isempty(fields)
                    for k=1:length(fields)     
                        if size(indin)==size(filters.(fields{k}))
                            indin=indin&filters.(fields{k});
                        else
                            disp('size of filters does not match. Call interfaces.LocData.filter before. Tried to fix it')
                            obj.filter(fields{k},layer)
                            indin=obj.inFilter(obj,layer,grouping);
                            filters
                            break       
                        end
                    end
                end
            end
        end
        
        function filter(obj,fields,layers,filtermode,minmax)
            %recaluclates all filters for specified layers.
            %filter(fields,layers,filtermode,filtervalues)
            %fields: cell array of fieldnames which to filter. If empty:
            %take all fields, also filter for channel
            %layers: all layers which to filter on. if empty: use all
            %layers
            %filtermode (optional): 'inlist' (filtervalues: list of integers to
            %keep),'minmax' (filtervalues: [minfilter maxfilter]). If not
            %specified: take from filter table
            
            if isempty(obj.loc)
                return
            end
            
            if nargin<3||isempty(layers) %filter all layers
                layers=1:obj.getPar('numberOfLayers');
            end
            
            if nargin<2||isempty(fields)%filter everything
                fields=fieldnames(obj.loc);
                 fields{end+1}='channel';
                 for k=1:length(layers)
                    obj.layer(layers(k)).groupfilter=[];
                    obj.layer(layers(k)).filter=[];
                 end
                 filterall=true;
            elseif ~iscell(fields)
                fields={fields};
                filterall=false;
            end

            for k=1:length(layers)  
                layerh=layers(k);
                s=obj.getPar('filtertable','layer',layerh);
                if ~isempty(s)
                fn=s(:,1);
                else
                    fn={};
                end
                for f=1:length(fields)
                    ind=find(strcmp(fn,fields{f}));
                    if (~isempty(ind) && (s{ind,7})) ||(~filterall&&isfield(obj.loc,fields{f})&&nargin>3) %if field specified: filter for sure
                        if nargin<4
                            minmax=s(ind,[2,6]);
                            filtermode='minmax';
                        end
                    elseif strcmp(fields{f},'channel')
                        filtermode='inlist';
                        if nargin <5
                            minmax=obj.getPar('channels','layer',k);
                        end
                    else 
                        continue
                    end   
                    obj.layer(layerh).filter.(fields{f})=anyfilter(obj.loc.(fields{f}),filtermode,minmax);
                    if ~isempty(obj.grouploc)
                        obj.layer(layerh).groupfilter.(fields{f})=anyfilter(obj.grouploc.(fields{f}),filtermode,minmax);
                    end
                end
            end            
        end
            
        function g=isgrouped(obj,layer)
            %returns grouping state of layer. Should belong to GuiChannel,
            %but link here is often useful
            %g=isgrouped(layer)
            g=zeros(length(layer),1);
            
            for k=1:length(layer)
                if layer(k)<=length(obj.layer)
             p=obj.getPar(['layer' num2str(layer(k)) '_']);
             g(k)=logical(p.groupcheck);  
                end
            end
            
            
        end
        
        function saveloc=savelocs(obj,filename,goodind,additionalsave)
            %saves localization data to specific filename. Only saves
            %ungroupd locs. Returns saved structure:
            %saveloc=savelocs(filename,goodind)
            %goodind (optional): indices which to save (applies to all
            %fields). Default: save all localizations.
            saveloc.loc=obj.loc;
            saveloc.file=obj.files.file;
            fieldsremove={'original_channel','groupindex','numberInGroup','colorfield'};
            saveloc.loc=myrmfield(saveloc.loc,fieldsremove);
            if nargin>2&&~isempty(goodind)
                fields=fieldnames(saveloc.loc);
                for k=1:length(fields)
                    saveloc.loc.(fields{k})=saveloc.loc.(fields{k})(goodind);
                end
            end
            if nargin>3&&~isempty(additionalsave)
                saveloc=copyfields(saveloc,additionalsave);
            end
            if nargin>1&&~isempty(filename)
                save(filename,'saveloc','-v7.3')
            end
        end
        
        function setLocData(obj,locData)
            %sets localizations to those from interfaces.LocalizationData object
            %setLocData(locData:interfaces.LocalizationData)
            obj.clear
            obj.loc=locData.loc;
            obj.files=locData.files;
            if isfield(locData,'SE')
                obj.SE=locData;
            end
        end
        
        function addLocData(obj,locData)  
            %adds localization from interfaces.LocalizationData object
            %addLocData(locData:interfaces.LocalizationData). Pad fields
            %which are not shared by both objects with zeros
            fn2=fieldnames(locData.loc);
            lennew=length(locData.loc.(fn2{1}));
            if ~isempty(obj.loc)
                fn1=fieldnames(obj.loc);
                lenold=length(obj.loc.(fn1{1}));
            else
                fn1={};
                lenold=0;
            end           
            %fill missing channels with 0
            fd=setdiff(fn2,fn1);
            for k=1:length(fd)
                obj.loc.(fd{k})=zeros(lenold,1,'like',locData.loc.(fd{k}));
            end  
            %add new data
            for k=1:length(fn2)
                obj.addloc(fn2{k},locData.loc.(fn2{k}))
            end     
            %add zeros for the missing data
            fd=setdiff(fn1,fn2);
            for k=1:length(fd)
                obj.addloc(fd{k},zeros(lennew,1,'like',obj.loc.(fd{k})));
            end                     
        end
        
        function sort(obj,varargin)
            %sorts localizations according to field1, field2,...
            %sort(field1,field2,...)
            sortm=zeros(length(obj.loc.xnm),length(varargin));
            for k=1:length(varargin)
                if isfield(obj.loc,varargin{k})
                    sortm(:,k)=(obj.loc.(varargin{k}));
                end
            end
            [~,sortind]=sortrows(sortm);
            fn=fieldnames(obj.loc);
            for k=1:length(fn)
                obj.loc.(fn{k})=obj.loc.(fn{k})(sortind);
            end
        end

        function locout= copy(obj,fields,indin,clearfilter)
            %creastes a deep copy of localization data.
            
            locout=interfaces.LocalizationData;
           locout.P=obj.P;
            locout.files=obj.files;
            
            if nargin<2
                locout.loc=obj.loc;
                locout.grouploc=obj.grouploc;

    %             locout.iscopy=true;
                for k=1:length(obj.layer)
                    locout.layer(k).filter=obj.layer(k).filter;
                    locout.layer(k).groupfilter=obj.layer(k).groupfilter;
                end
            else
                if nargin<4
                    clearfilter=false;
                end
                    
                fn=fieldnames(obj.loc);
                if isempty(fn)
                    return;
                end
                
                numlocs=length(obj.loc.(fn{1}));
                numglocs=length(obj.grouploc.(fn{1}));
                
                if nargin<3
                    indu=true(numlocs,1);
                    indg=true(numglocs,1);
                else
                    if length(indin)==numlocs;
                        indu=indin;
                        indg=indungrouped2grouped(indin,obj.loc.groupindex);
                    else
                        indg=indin;
                        indu=indgrouped2ungrouped(indin,obj.loc.groupindex);
                    end
                end
                if isempty(fields)
                    fields=fn;
                else
                    fields=intersect(fields,fn);
                end
                for k=1:length(fields)
                    locout.loc.(fields{k})=obj.loc.(fields{k})(indu);
                end
                for k=1:length(fields)
                    locout.grouploc.(fields{k})=obj.grouploc.(fields{k})(indg);
                end
%                 fields=fieldnames(obj.grouploc);
                if clearfilter
                    for l=1:length(obj.layer);
                    locout.layer(l).filter=[];
                    locout.layer(l).groupfilter=[];
                    end
                else
                    

    %                 l=length(obj.layer);
                    for l=1:length(obj.layer);
                        fields=fieldnames(obj.layer(l).filter);
                        for k=1:length(fields)
                            locout.layer(l).filter.(fields{k})=obj.layer(l).filter.(fields{k})(indu);
                        end        
                        fields=fieldnames(obj.layer(l).groupfilter);
                        for k=1:length(fields)
                            locout.layer(l).groupfilter.(fields{k})=obj.layer(l).groupfilter.(fields{k})(indg);
                        end                
                    end
                end
            end
        end
        
        function regroup(obj,dx,dt)  
            %groups localization data.  
            %regroup(dx,dt): dx: maximum distance to be grouped. dt:
            %maximum number of frames localzaiton can be dark. if not
            %specified: taken from global parameters 'group_dx' and
            %'group_dt'
            if nargin<2
                dx=obj.getPar('group_dx');
            end
            if nargin<3
                dt=obj.getPar('group_dt');
            end
            
%             obj.setPar('status','grouping')
            grouper=Grouper;
            grouper.attachLocData(obj);
            grouper.connect(dx,dt,'frame','xnm','ynm','filenumber','channel')
            grouper.combine;  
            obj.filter;
%             locfields=fieldnames(obj.loc);
%             if ~obj.iscopy
%                 obj.setPar('locFields',locfields);
%             end
%             obj.setPar('status','grouping done')
        end    
    end    
end