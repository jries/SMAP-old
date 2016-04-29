function initGuiAfterLoad(obj)
if isempty(obj.locData.loc)
    return
end
obj.locData.regroup;
% obj.locData.filter;
% fmax=max(obj.locData.loc.frame);
%this calls histogram maker etc
% obj.setPar('frame_max',fmax);
% obj.setPar('frame_min',1);

fl={obj.locData.files.file(:).name};
for k=1:length(fl)
    flh=strrep(fl{k},'\',filesep);
    flh=strrep(flh,'/',filesep);
    [~,fls{k}]=fileparts(flh);
end

obj.setPar('filelist_long',fl,'String');
% obj.locData.filter;
obj.setPar('cam_pixelsize_nm',obj.locData.files.file(1).info.pixsize*1000)
% obj.setPar('mainfile',fl{1});

locfields=fieldnames(obj.locData.loc);
obj.setPar('locFields',locfields,'String');

obj.setPar('filelist_short',fls,'String');
obj.locData.filter;
obj.setPar('currentfileinfo',obj.locData.files.file(1).info)



if ~isempty(strfind(obj.getPar('mainfile'),'_sml'))
   tg=obj.getPar('mainGui').guihandles.maintab;
   tg.SelectedTab=tg.Children(3);
end
end