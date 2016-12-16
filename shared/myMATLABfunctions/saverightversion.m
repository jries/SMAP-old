function v=saverightversion(file,ls)
v='v7';
save(file,'-struct','ls','-v7');
[msg,msgid]=lastwarn;
if strcmp(msgid,'MATLAB:save:sizeTooBigForMATFile')
    save(file,'-struct','ls','-v7.3');
    disp('File is now being saved as v7.3');
    v='v7.3';
    lastwarn('cleared');
end


