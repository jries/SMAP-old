function sout=loadmatparts(f)

% [path,file,ext]=fileparts(f);
ind=regexp(f,'_p\d*.mat');
if ~isempty(ind)
f=[f(1:ind-1) '.mat'];
end

ld=load(f);
if ~isfield(ld,'lds') && ~isfield(ld,'S') %not part file
    sout=ld;
    return
end

lh=subsref(ld.lds,ld.S);
for k=1:length(ld.partnames)
    lk=load(ld.partnames{k});
    lh=copyfields(lh,lk.ltemp);
end
sout=ld.lds;
sout=subsasgn(sout,ld.S,lh);
end