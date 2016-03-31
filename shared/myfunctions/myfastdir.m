function outs=myfastdir(pfad,search)
flist=evalc('dir(pfad)');
s2=strrep(search,'*','[a-z_0-9]+');
s2=strrep(s2,'.','\.');
[ind ind2]=regexpi(flist,s2,'start','end');
out{length(ind)}=flist(ind(1):ind2(1));
for k=1:length(ind)
    out{k}=flist(ind(k):ind2(k));
end

outs=sort(out);