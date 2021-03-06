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
[path,file]=fileparts(f);

lh=subsref(ld.lds,ld.S);
for k=1:length(ld.partnames)
    [~,fileh,ext]=fileparts(ld.partnames{k});
    lk=load([path filesep fileh ext]);
    lh=copyfields(lh,lk.ltemp);
end
sout=ld.lds;
sout=subsasgn(sout,ld.S,lh);
end