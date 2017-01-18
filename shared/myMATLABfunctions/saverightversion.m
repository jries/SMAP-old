function v=saverightversion(file,ls,version)
if nargin <3 || isempty(version)
    v='-v7';
else
    v=version;
end
save(file,'-struct','ls',v);
[msg,msgid]=lastwarn;
if strcmp(msgid,'MATLAB:save:sizeTooBigForMATFile')
    save(file,'-struct','ls','-v7.3');
    disp('File is now being saved as v7.3');
    v='-v7.3';
    lastwarn('cleared');
end


