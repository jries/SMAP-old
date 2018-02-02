% removeFieldsSML

fieldstosave={'phot','frame','PSFxnm','locprecnm','xnm','ynm','bg','channel','filenumber','dummy'};

startdirectory='c:/Users/*_sml.mat';
h=selectManyFiles(startdirectory);
waitfor(h.handle)

filelist=h.filelist;
for k=1:length(filelist)
    l=load(filelist{k});
    lout=l;
    fn=fieldnames(l.saveloc.loc);
    fieldsremove=setdiff(fn,fieldstosave);
    notfound=setdiff(fieldstosave,fn);
    if ~isempty(notfound)
        disp(['fields not found: ' notfound{1}]);
    end
    lout.saveloc.loc=rmfield(l.saveloc.loc,fieldsremove);
    newname=strrep(filelist{k},'_sml.mat','_s_sml.mat');
    v=saverightversion(file,lout,'-v7');
    disp(v)
end