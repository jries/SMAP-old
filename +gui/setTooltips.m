function setTooltips(obj)
%GuiFile
o=obj.children.guiFile.guihandles;
o.load.TooltipString='load localization data or image. It is recommended to load at least one localization data before an image.';
o.add.TooltipString='add loclaization data or image';
o.remove.TooltipString='remove file';
o.savemode.TooltipString='select what to save';
o.group_b.TooltipString='group localizations in consecutive frames';
o.group_dx.TooltipString=sprintf('distance in nm which two locs can be apart \n and still grouped together');
o.group_dt.TooltipString=sprintf('number of frames locs can be missing \n and still grouped together');

%status
% o=obj.children.status.guihandles;
% o.profile.TooltipString='set profiler on or off. Standard should be off';
% o.profileb.TooltipString='show profile results and clear profile history';

% layers
% layernames={'l1', 'l2','l3','l4','l5','l6'};
% for k=1:length(layernames)
% o=obj.children.guiRender.children.(layernames{k}).guihandles;
% pard.layercheck.TooltipString='switch layer on and off';
% pard.channels.TooltipString='channels to display. Use a,b,c and a:c notation';
% pard.rendermode.TooltipString='how to render image. DL is diffraction limited';
% pard.renderp1.TooltipString=sprintf('normal: intensity coded, z: z color coded, \n param: select field which to color code');
% pard.renderfield.TooltipString='field to color code';
% pard.groupcheck.TooltipString='use grouped or ungrouped los';
% pard.ch_filelist.TooltipString='which file (loc or image) to display';
% pard.imaxtoggle.TooltipString='toggle absolute intensity maximum (Imax) or quantile';
% pard.imax.TooltipString='absolut intensity or quantile (0<q<1) or v for q=1-10^(v), v<0';
% pard.cb.TooltipString='range of values to fill the lookup table';
% pard.colorfield_min.TooltipString=pard.cb.TooltipString;
% pard.colorfield_max.TooltipString=pard.cb.TooltipString;
% pard.remout.TooltipString=sprintf('if checked: remove loclizations outside lut. \n If unchecked: set them to maximum color');
% pard.lut.TooltipString='select the lookup table';
% pard.parbutton.TooltipString='Additional render paramters';
% 
% hist=obj.children.guiRender.children.(layernames{k}).children.histgui.guihandles;
% hist.lockrange.TooltipString='if checked, the range will be restricted to the width given here';
% hist.range.TooltipString='width of range of values';
% hist.showall.TooltipString='if checked, this filter has no effect and all values are used';
% 
% end 

%format
o=obj.children.guiRender.children.guiFormat.guihandles;
% o.hlayers.TooltipString='select which layers to display';
o.hplus.TooltipString='zoom out, increase pixelsize';
o.hminus.TooltipString='zoom in, decrease pixelsize';
o.pixrec.TooltipString='pixel size for reconstruction (nm)';
o.pref1.TooltipString=sprintf('preset pixel size. To change this value:  \n 1. click on button and leave mouse over it \n 2. type pixelsize in nm, finishe with Enter');
o.pref2.TooltipString=o.pref1.TooltipString;
o.resetview.TooltipString='adjust pixelsize to fit all width or height';
o.parformat.TooltipString='additional global parameters for rendering';
o.redrawov.TooltipString='redraw overview image using the settings of the layer tab: ovim 6';
o.overview_select.TooltipString='toggels between overview image and histogram view';
o.linewidth_roi.TooltipString='width of the roi when using the line';
o.roi4.TooltipString='line roi, width set below';
o.roi5.TooltipString='point roi, rectangular roi around with width and height set below';
o.roi1.TooltipString='rectangular roi';
o.roi2.TooltipString='circular roi';
o.roi6.TooltipString='polynome roi';
o.roi3.TooltipString='free roi';
