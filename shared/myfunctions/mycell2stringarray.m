function sout=mycell2stringarray(cin)
len=max(cellfun(@length,cin));
sout(length(cin),len)='a';
for k=1:length(cin)
    sout(k,1:length(cin{k}))=cin{k};
end
end