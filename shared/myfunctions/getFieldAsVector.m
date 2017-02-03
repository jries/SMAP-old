function v=getFieldAsVector(p,varargin)
% ls=length(p);
isarray=false;
for k=length(p):-1:1
    try
    if iscell(p)
%         if isfield(p{k},varargin{1})
            vh=p{k}.(varargin{1});
%         end
%         isarray=true;
    else
%         if isfield(p(k),varargin{1})
            vh=p(k).(varargin{1});
%         end
    end
    for f=2:length(varargin)
%         if ~isempty(vh)
            nv=vh.(varargin{f});
            vh=nv;
%         end
    end
    if isarray||(numel(vh)==1 && (isnumeric(vh)||islogical(vh)))
        v(k)=vh;
        isarray=true;
    else
        v{k}=vh;
    end
    catch err
    if isarray
        v(k)=NaN;
    else
        v{k}=NaN;
    end
    end

end