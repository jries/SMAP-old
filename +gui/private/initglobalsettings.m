function initglobalsettings(obj)
%defines global settings

obj.createGlobalSetting('guiPluginConfigFile','GUI','Configuration file for GUI plugin structure. Delete path and save to reset plugins.',struct('Style','file','String','settings/SimpleGui.txt'))
obj.createGlobalSetting('customMenuFile','GUI','Configuration file for custom menu. Delete path and save to not have any custom menu.',struct('Style','file','String',' '))
obj.createGlobalSetting('mainLocalizeWFFile','GUI','Description file for fitting workflow, e.g. settings/workflows/fit_tif_wavelet.txt',struct('Style','file','String','settings/workflows/fit_any_wavelet.txt'))
object=struct('Style','saveparameter','String','save current gui plugin configuration');
 obj.createGlobalSetting('guimodules','GUI','',object);
 
if ispc
    fijipath='C:/Program Files/Fiji/scripts';
else
    fijipath='/Applications/Fiji.app/scripts';
end
obj.createGlobalSetting('MMpath','Directories','The directory of Micro-Manager in which ij.jar is found:',struct('Style','dir','String','MMpath'))   
obj.createGlobalSetting('fijipath','Directories','The directory of /Fiji/scripts:',struct('Style','dir','String',fijipath))

structure=struct('Style','checkbox','String','Ask for Auto Check on','Value',1);
obj.createGlobalSetting('SE_autosavecheck','ROIManager','Ask for autosave on when starting ROI manager',structure);
 obj.createGlobalSetting('saveas73','File','.mat format. -v7.3: slow and larger files, but ompatibel with >2GB',...
    struct('Style','checkbox','String','Always use -v7.3','Value',0));