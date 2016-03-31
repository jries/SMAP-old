function v=getFieldAsVector(p,varargin)
% ls=length(p);
isarray=false;
for k=length(p):-1:1
    if iscell(p)
        vh=p{k}.(varargin{1});
%         isarray=true;
    else
        vh=p(k).(varargin{1});
    end
    for f=2:length(varargin)
        nv=vh.(varargin{f});
        vh=nv;
    end
    if isarray||(numel(vh)==1 && isnumeric(vh))
        v(k)=vh;
        isarray=true;
    else
        v{k}=vh;
    end

end