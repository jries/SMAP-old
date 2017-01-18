function mij=openfiji(obj)



mij=obj.getPar('MIJ');
if isempty(mij) %open fiji
    
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