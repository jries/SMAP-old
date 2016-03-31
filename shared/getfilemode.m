
function [mode, emptylocs]=getfilemode(file)
[~,f,ext]=fileparts(file);
mode='';
emptylocs=true;
if strcmpi(ext,'.tif')
    mode='tif';
    emptylocs=false;
    return
end
if strcmpi(ext,'.mat')&& (strcmp(f(end-4:end),'_sml'))
    mode='sml';
    emptylocs=true;
    return
end
if strcmpi(ext,'.csv')
    mode='csv';
    emptylocs=true;
    return
end
if strcmpi(ext,'.mat')&& ~isempty(strfind(f,'_sites'))
    mode='sites';
    emptylocs=true;
    return
end
if ~isempty(strfind(f,'fitpos'))
    mode='fitpos';
    emptylocs=true;
    return
end
vars=whos('-file',file);
if sum(strcmpi({vars.name},'fileformat'))
    f=load(file,'fileformat');
    mode=f.fileformat.name;
    if any(strcmp(mode,{'guiparameters'}))
        emptylocs=false;
    end
    return
end

if sum(strcmpi({vars.name},'fitpos'))
    mode='fitpos';
    return
end

if sum(strcmpi({vars.name},'saveloc'))
    mode='sml';
    return
end

if sum(strcmpi({vars.name},'SEsettings'))
    mode='settings';
    return
end

end
