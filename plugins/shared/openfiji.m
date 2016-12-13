function mij=openfiji(obj)



mij=obj.getPar('MIJ');
if isempty(mij) %open fiji
    if ispc
        fijipath='C:/Program Files/Fiji/scripts';
    else
        fijipath='/Applications/Fiji.app/scripts';
    end
    obj.createGlobalSetting('fijipath','Directories2','The directory of /Fiji/scripts:',struct('Style','dir','String',fijipath))
    fijipath=obj.getGlobalSetting('fijipath');  
    if ~exist(fijipath,'dir')
        mij=[];
        errordlg('cannot find Fiji, please select Fiji directory in menu SMAP/Preferences...')
        return
    end
    
    dir=pwd;
    obj.setPar('status','open Fiji');
    addpath(fijipath)
    Miji();
    mij=MIJ;
    cd(dir);
    obj.setPar('MIJ',mij);
end
end