classdef Grouper< interfaces.LocDataInterface
    properties
        combinemodes
        indsortlist
    end
    methods
        function obj=Grouper(varargin)            
            if nargin>0
                obj.attachLocData(varargin{1});
            end
%             if nargin>1
%                 obj.attachPar(varargin{2});
%             end
%             obj.inputParameters={'group_dx','group_dt'};
            obj.combinemodes.frame='first';
            obj.combinemodes.xnm='mean';
            obj.combinemodes.ynm='mean';
            obj.combinemodes.znm='mean';
            obj.combinemodes.bg='sum';
            obj.combinemodes.phot='sum';
            obj.combinemodes.PSFxnm='square';
            obj.combinemodes.PSFynm='square';
            obj.combinemodes.locprecnm='locp';
            obj.combinemodes.locprecxnm='locp';
            obj.combinemodes.locprecynm='locp';
            obj.combinemodes.locprecznm='locp';
            obj.combinemodes.CRLBphot='min';
            obj.combinemodes.frame='first';
            obj.combinemodes.channel='first';
            obj.combinemodes.logLikelihood='max';
            obj.combinemodes.loglikelihood='max';
            obj.combinemodes.numberInGroup='first';
            obj.combinemodes.groupindex='first';
            obj.combinemodes.filenumber='first';
            obj.combinemodes.neighbours='sum';
        end
        function connect(obj,dx,dt,framef,xf,yf,varargin)
%             obj.status('group localizations')
            %here not clear. like this or exchanged?
            x=double(obj.locData.getloc(xf).(xf));
            y=double(obj.locData.getloc(yf).(yf));
            frame=double(obj.locData.getloc(framef).(framef));
            
            sortmatrix=([frame,x,y]);  %modify:unconnected
            for k=1:length(varargin)
                sortmatrix=[double(obj.locData.getloc(varargin{k}).(varargin{k})),sortmatrix];
            end
            sm=size(sortmatrix);
            [~,indsort]=sortrows(sortmatrix,1:sm(2));
            
            maxactive=10000;
            list=connectsingle2c(double(x(indsort)),double(y(indsort)),double(frame(indsort)),double(dx),int32(dt),int32(maxactive));
%             listm=connectsingle2mat(double(x(indsort)),double(y(indsort)),double(frame(indsort)),double(dx),int32(dt),int32(maxactive));
            numbers=1:sm(1);
            indold=numbers(indsort);
            [~,indback]=sort(indold);
            listback=list(indback);
            obj.locData.setloc('groupindex',listback)
            
             %number of locs
            [listsort,indsort2]=sort(list);
%             numbergroup=zeros(size(listsort));
%             whos listsort
            if ~isempty(listsort)
            numbergroup=countlocs(double(listsort));
            
           
%             groupc=listsort(1);
%             ng=0;
%             inds=1;
%             for k=1:length(numbergroup)
%                 if listsort(k)~=groupc
%                     groupc=listsort(k);
%                     numbergroup(inds:k-1)=ng;
%                     ng=0;
%                     inds=k;
%                 end
%                 ng=ng+1;
% %                 numbergroup(k)=ng;
%             end
            
            
            indold2=indold(indsort2);
            [~,indback2]=sort(indold2);
            obj.locData.setloc('numberInGroup',numbergroup(indback2));
            [~,obj.indsortlist]=sort(listback);
            end
%             obj.status('group localizations done')

            
        end
        function combine(obj,field,combinemode,weights) %no field etc: group everythign for which we have combinemodes
            
%             obj.status('group: combine fields')
            if nargin==1 %do all
%                 fn=fieldnames(obj.combinemodes);
                fn2=fieldnames(obj.locData.loc);
                if isempty(obj.locData.loc.(fn2{1}))  
                   return
                end
%                 fnall=intersect(fn,fn2);
                weights=1./(obj.locData.getloc('locprecnm').locprecnm);
                if isempty(weights)
                    weights=ones(size(obj.locData.loc.(fn2{1})));  
                else
                    weights(isinf(weights))=1;
                end
                for k=1:length(fn2)
                    if isfield(obj.combinemodes,fn2{k})
                        combinemode=obj.combinemodes.(fn2{k});
                    else
                        combinemode='mean';
                    end
                    obj.combine(fn2{k},combinemode,weights);
                end
                
            else
                
            if nargin>2
                obj.combinemodes.(field)=combinemode;
            elseif  nargin==2
                if isfield(obj.combinemodes,'field')
                combinemode=obj.combinemodes.(field);
                else
                    combinemode='mean';
                end
            end
            
            v=obj.locData.getloc(field).(field);
          
            list=obj.locData.getloc('groupindex').groupindex;
            
            if nargin <4
                weights=ones(size(v));              
            end          
            mode=0;
            switch combinemode
                case 'sum'
                    v2=v;  
                case 'mean'
                    v2=v.*weights;                 
                case 'square'
                    v2=v.^2.*weights.^2;
                    weights=weights.^2;
                case 'min'
                    v2=-v;
                    mode=1;
                case 'max'
                    v2=v;
                    mode=1;
                case 'first'
                    mode=2;
                    v2=v;   
                case 'locp'
                    v2=1./v.^2;  
            end
            
            

            %combine modes
            %<x>
            %sqrt(<x^2>)
            %first
            %sum
%             [listsort,indsort]=sort(list);
            indsort=obj.indsortlist;
            listsort=list(indsort);
            vwsort=v2(indsort);

            if mode==2
                dl=[diff(listsort);1];
                vwout2=vwsort(dl==1);
%                 listsort(end)
%                 sum(dl==1)
            else 
                if mode==0
                vwout=sumcombine(double(vwsort),double(listsort));
                gweights=sumcombine(double(weights(indsort)),double(listsort));
                switch combinemode
                    case 'sum'
                        vwout2=vwout;  
                    case 'mean'
                        vwout2=vwout./gweights;                 
                    case 'square'
                        vwout2=sqrt(vwout./gweights); 
                    case 'locp'
                        vwout2=1./sqrt(vwout);
                end   
                elseif mode==1
                    vwout=maxcombine(double(vwsort),double(listsort));
                    if strcmp(combinemode,'min')
                        vwout2=-vwout;
                    else
                        vwout2=vwout;
                    end
                end
            end
           
            obj.locData.grouploc.(field)=cast(vwout2,'like',v);   
            
            end
%             obj.status('group: combine fields done')
        end
    end
end
