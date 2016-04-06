function out=plugins(fn1,fn2,fn3,varargin) 
if nargin>0 
fstr=[fn1 '.' fn2 '.' fn3]; 
switch fstr 
case 'Analyze.counting.cfithist' 
   module=Analyze.counting.cfithist(varargin{:}); 
case 'Analyze.counting.cloadhist' 
   module=Analyze.counting.cloadhist(varargin{:}); 
case 'Analyze.counting.cmakehist' 
   module=Analyze.counting.cmakehist(varargin{:}); 
case 'Analyze.counting.csegmentcluster' 
   module=Analyze.counting.csegmentcluster(varargin{:}); 
case 'Analyze.counting.locsfromSE' 
   module=Analyze.counting.locsfromSE(varargin{:}); 
case 'Analyze.measure.FRCresolution' 
   module=Analyze.measure.FRCresolution(varargin{:}); 
case 'Analyze.measure.Locstatistics' 
   module=Analyze.measure.Locstatistics(varargin{:}); 
case 'Analyze.measure.PSFanalysis' 
   module=Analyze.measure.PSFanalysis(varargin{:}); 
case 'Analyze.measure.StatsVsTime' 
   module=Analyze.measure.StatsVsTime(varargin{:}); 
case 'Analyze.measure.clusterMIiSR' 
   module=Analyze.measure.clusterMIiSR(varargin{:}); 
case 'Analyze.measure.lineprofile' 
   module=Analyze.measure.lineprofile(varargin{:}); 
case 'Analyze.other.BlinkingMoviePresentation' 
   module=Analyze.other.BlinkingMoviePresentation(varargin{:}); 
case 'Analyze.other.CompareToGroundTruth' 
   module=Analyze.other.CompareToGroundTruth(varargin{:}); 
case 'Analyze.other.PushToFiji' 
   module=Analyze.other.PushToFiji(varargin{:}); 
case 'Analyze.other.ShowHistory' 
   module=Analyze.other.ShowHistory(varargin{:}); 
case 'Analyze.other.VersatileRenderer' 
   module=Analyze.other.VersatileRenderer(varargin{:}); 
case 'Analyze.other.density_calculator' 
   module=Analyze.other.density_calculator(varargin{:}); 
case 'Analyze.projects.Platelets_Band' 
   module=Analyze.projects.Platelets_Band(varargin{:}); 
case 'Analyze.render.VoronoiRenderer' 
   module=Analyze.render.VoronoiRenderer(varargin{:}); 
case 'Analyze.sr3D.CheckCalibration' 
   module=Analyze.sr3D.CheckCalibration(varargin{:}); 
case 'Analyze.sr3D.Sxsy2z' 
   module=Analyze.sr3D.Sxsy2z(varargin{:}); 
case 'Analyze.sr3D.VISPexport' 
   module=Analyze.sr3D.VISPexport(varargin{:}); 
case 'Analyze.sr3D.Viewer3DV01' 
   module=Analyze.sr3D.Viewer3DV01(varargin{:}); 
case 'Analyze.sr3D.calibrateAstig' 
   module=Analyze.sr3D.calibrateAstig(varargin{:}); 
case 'Analyze.sr3D.calibrateAstigDeep' 
   module=Analyze.sr3D.calibrateAstigDeep(varargin{:}); 
case 'Analyze.sr3D.sideview' 
   module=Analyze.sr3D.sideview(varargin{:}); 
case 'Analyze.sr3D.volume3D' 
   module=Analyze.sr3D.volume3D(varargin{:}); 
case 'File.Load.Loader_auto' 
   module=File.Load.Loader_auto(varargin{:}); 
case 'File.Load.Loader_csv' 
   module=File.Load.Loader_csv(varargin{:}); 
case 'File.Load.Loader_ioPlugin' 
   module=File.Load.Loader_ioPlugin(varargin{:}); 
case 'File.Load.Loader_settings' 
   module=File.Load.Loader_settings(varargin{:}); 
case 'File.Load.Loader_sml' 
   module=File.Load.Loader_sml(varargin{:}); 
case 'File.Load.Loader_tif' 
   module=File.Load.Loader_tif(varargin{:}); 
case 'File.Load.Loader_workflow' 
   module=File.Load.Loader_workflow(varargin{:}); 
case 'File.Save.SMLMsaver' 
   module=File.Save.SMLMsaver(varargin{:}); 
case 'File.Save.TifSaver' 
   module=File.Save.TifSaver(varargin{:}); 
case 'File.Save.saver_csv' 
   module=File.Save.saver_csv(varargin{:}); 
case 'File.Save.saver_settings' 
   module=File.Save.saver_settings(varargin{:}); 
case 'Process.Assign2C.Get2CIntImagesWF' 
   module=Process.Assign2C.Get2CIntImagesWF(varargin{:}); 
case 'Process.Assign2C.Get2CIntLoc' 
   module=Process.Assign2C.Get2CIntLoc(varargin{:}); 
case 'Process.Assign2C.Intensity2Channel' 
   module=Process.Assign2C.Intensity2Channel(varargin{:}); 
case 'Process.Assign2C.SALMInt2z' 
   module=Process.Assign2C.SALMInt2z(varargin{:}); 
case 'Process.Drift.driftcorrection' 
   module=Process.Drift.driftcorrection(varargin{:}); 
case 'Process.Drift.showdriftinfo' 
   module=Process.Drift.showdriftinfo(varargin{:}); 
case 'Process.Images.Tiff_remove_BG' 
   module=Process.Images.Tiff_remove_BG(varargin{:}); 
case 'Process.Modify.Connect2Unconnect' 
   module=Process.Modify.Connect2Unconnect(varargin{:}); 
case 'Process.Modify.MathParser' 
   module=Process.Modify.MathParser(varargin{:}); 
case 'Process.Modify.Mathematics' 
   module=Process.Modify.Mathematics(varargin{:}); 
case 'Process.Modify.Mirror' 
   module=Process.Modify.Mirror(varargin{:}); 
case 'Process.Modify.RemoveLocs' 
   module=Process.Modify.RemoveLocs(varargin{:}); 
case 'Process.Modify.RenameFields' 
   module=Process.Modify.RenameFields(varargin{:}); 
case 'Process.Register.ApplyTransform' 
   module=Process.Register.ApplyTransform(varargin{:}); 
case 'Process.Register.RegisterImages' 
   module=Process.Register.RegisterImages(varargin{:}); 
case 'Process.Register.RegisterLocs2' 
   module=Process.Register.RegisterLocs2(varargin{:}); 
case 'ROIManager.Analyze.AnalyzeCME2Cside' 
   module=ROIManager.Analyze.AnalyzeCME2Cside(varargin{:}); 
case 'ROIManager.Analyze.AnalyzeRingCME' 
   module=ROIManager.Analyze.AnalyzeRingCME(varargin{:}); 
case 'ROIManager.Analyze.PlotRingCME' 
   module=ROIManager.Analyze.PlotRingCME(varargin{:}); 
case 'ROIManager.Analyze.RemoveEmptyCells' 
   module=ROIManager.Analyze.RemoveEmptyCells(varargin{:}); 
case 'ROIManager.Analyze.export_tiffs' 
   module=ROIManager.Analyze.export_tiffs(varargin{:}); 
case 'ROIManager.Evaluate.CME2CSide' 
   module=ROIManager.Evaluate.CME2CSide(varargin{:}); 
case 'ROIManager.Evaluate.CME2DRing' 
   module=ROIManager.Evaluate.CME2DRing(varargin{:}); 
case 'ROIManager.Evaluate.CME2DRingFit' 
   module=ROIManager.Evaluate.CME2DRingFit(varargin{:}); 
case 'ROIManager.Evaluate.CME3DDSpherefit' 
   module=ROIManager.Evaluate.CME3DDSpherefit(varargin{:}); 
case 'ROIManager.Evaluate.CME3DDSpherefit_DC' 
   module=ROIManager.Evaluate.CME3DDSpherefit_DC(varargin{:}); 
case 'ROIManager.Evaluate.countingStatistics' 
   module=ROIManager.Evaluate.countingStatistics(varargin{:}); 
case 'ROIManager.Evaluate.generalStatistics' 
   module=ROIManager.Evaluate.generalStatistics(varargin{:}); 
case 'ROIManager.Segment.segmentCellsCME' 
   module=ROIManager.Segment.segmentCellsCME(varargin{:}); 
case 'ROIManager.Segment.segmentNPC' 
   module=ROIManager.Segment.segmentNPC(varargin{:}); 
case 'ROIManager.Segment.segmentSitesPolarCME' 
   module=ROIManager.Segment.segmentSitesPolarCME(varargin{:}); 
case 'WorkflowModules.Filters.BG_wavelet' 
   module=WorkflowModules.Filters.BG_wavelet(varargin{:}); 
case 'WorkflowModules.Filters.DisplayChooser' 
   module=WorkflowModules.Filters.DisplayChooser(varargin{:}); 
case 'WorkflowModules.Filters.ImageFilter' 
   module=WorkflowModules.Filters.ImageFilter(varargin{:}); 
case 'WorkflowModules.Filters.ImageNormalize' 
   module=WorkflowModules.Filters.ImageNormalize(varargin{:}); 
case 'WorkflowModules.Filters.MedianBGcalculator' 
   module=WorkflowModules.Filters.MedianBGcalculator(varargin{:}); 
case 'WorkflowModules.Filters.fitOnBackground' 
   module=WorkflowModules.Filters.fitOnBackground(varargin{:}); 
case 'WorkflowModules.Fitters.EMCCD_SE_MLE_GPU' 
   module=WorkflowModules.Fitters.EMCCD_SE_MLE_GPU(varargin{:}); 
case 'WorkflowModules.Fitters.FalconHD2D' 
   module=WorkflowModules.Fitters.FalconHD2D(varargin{:}); 
case 'WorkflowModules.Fitters.RadialSymmetry2D' 
   module=WorkflowModules.Fitters.RadialSymmetry2D(varargin{:}); 
case 'WorkflowModules.Fitters.RadialSymmetry3D' 
   module=WorkflowModules.Fitters.RadialSymmetry3D(varargin{:}); 
case 'WorkflowModules.Fitters.fitterGUI' 
   module=WorkflowModules.Fitters.fitterGUI(varargin{:}); 
case 'WorkflowModules.IntensityCalculator.EvaluateIntensity' 
   module=WorkflowModules.IntensityCalculator.EvaluateIntensity(varargin{:}); 
case 'WorkflowModules.IntensityCalculator.IntLoc2pos' 
   module=WorkflowModules.IntensityCalculator.IntLoc2pos(varargin{:}); 
case 'WorkflowModules.IntensityCalculator.pushLocs' 
   module=WorkflowModules.IntensityCalculator.pushLocs(varargin{:}); 
case 'WorkflowModules.IntensityCalculator.roi2int_fitG' 
   module=WorkflowModules.IntensityCalculator.roi2int_fitG(varargin{:}); 
case 'WorkflowModules.IntensityCalculator.roi2int_sumG' 
   module=WorkflowModules.IntensityCalculator.roi2int_sumG(varargin{:}); 
case 'WorkflowModules.Loaders.CameraConverter' 
   module=WorkflowModules.Loaders.CameraConverter(varargin{:}); 
case 'WorkflowModules.Loaders.GrabFijiStacks' 
   module=WorkflowModules.Loaders.GrabFijiStacks(varargin{:}); 
case 'WorkflowModules.Loaders.TifLoader' 
   module=WorkflowModules.Loaders.TifLoader(varargin{:}); 
case 'WorkflowModules.Loaders.workflowstarter' 
   module=WorkflowModules.Loaders.workflowstarter(varargin{:}); 
case 'WorkflowModules.LocProcessors.LocFilter' 
   module=WorkflowModules.LocProcessors.LocFilter(varargin{:}); 
case 'WorkflowModules.LocProcessors.LocSaver' 
   module=WorkflowModules.LocProcessors.LocSaver(varargin{:}); 
case 'WorkflowModules.LocProcessors.OnlineReconstruction' 
   module=WorkflowModules.LocProcessors.OnlineReconstruction(varargin{:}); 
case 'WorkflowModules.LocProcessors.PlotLocsPreview' 
   module=WorkflowModules.LocProcessors.PlotLocsPreview(varargin{:}); 
case 'WorkflowModules.Peakfinders.PeakFinder' 
   module=WorkflowModules.Peakfinders.PeakFinder(varargin{:}); 
case 'WorkflowModules.Peakfinders.PeakFinderNMSWF' 
   module=WorkflowModules.Peakfinders.PeakFinderNMSWF(varargin{:}); 
case 'WorkflowModules.Peakfinders.RoiAdder' 
   module=WorkflowModules.Peakfinders.RoiAdder(varargin{:}); 
case 'WorkflowModules.Peakfinders.RoiCutterWF' 
   module=WorkflowModules.Peakfinders.RoiCutterWF(varargin{:}); 
end 
module.pluginpath={fn1,fn2,fn3}; 
out=module; 
else 
out.Analyze.counting.cfithist={'Analyze','counting','cfithist'}; 
out.Analyze.counting.cloadhist={'Analyze','counting','cloadhist'}; 
out.Analyze.counting.cmakehist={'Analyze','counting','cmakehist'}; 
out.Analyze.counting.csegmentcluster={'Analyze','counting','csegmentcluster'}; 
out.Analyze.counting.locsfromSE={'Analyze','counting','locsfromSE'}; 
out.Analyze.measure.FRCresolution={'Analyze','measure','FRCresolution'}; 
out.Analyze.measure.Locstatistics={'Analyze','measure','Locstatistics'}; 
out.Analyze.measure.PSFanalysis={'Analyze','measure','PSFanalysis'}; 
out.Analyze.measure.StatsVsTime={'Analyze','measure','StatsVsTime'}; 
out.Analyze.measure.clusterMIiSR={'Analyze','measure','clusterMIiSR'}; 
out.Analyze.measure.lineprofile={'Analyze','measure','lineprofile'}; 
out.Analyze.other.BlinkingMoviePresentation={'Analyze','other','BlinkingMoviePresentation'}; 
out.Analyze.other.CompareToGroundTruth={'Analyze','other','CompareToGroundTruth'}; 
out.Analyze.other.PushToFiji={'Analyze','other','PushToFiji'}; 
out.Analyze.other.ShowHistory={'Analyze','other','ShowHistory'}; 
out.Analyze.other.VersatileRenderer={'Analyze','other','VersatileRenderer'}; 
out.Analyze.other.density_calculator={'Analyze','other','density_calculator'}; 
out.Analyze.projects.Platelets_Band={'Analyze','projects','Platelets_Band'}; 
out.Analyze.render.VoronoiRenderer={'Analyze','render','VoronoiRenderer'}; 
out.Analyze.sr3D.CheckCalibration={'Analyze','sr3D','CheckCalibration'}; 
out.Analyze.sr3D.Sxsy2z={'Analyze','sr3D','Sxsy2z'}; 
out.Analyze.sr3D.VISPexport={'Analyze','sr3D','VISPexport'}; 
out.Analyze.sr3D.Viewer3DV01={'Analyze','sr3D','Viewer3DV01'}; 
out.Analyze.sr3D.calibrateAstig={'Analyze','sr3D','calibrateAstig'}; 
out.Analyze.sr3D.calibrateAstigDeep={'Analyze','sr3D','calibrateAstigDeep'}; 
out.Analyze.sr3D.sideview={'Analyze','sr3D','sideview'}; 
out.Analyze.sr3D.volume3D={'Analyze','sr3D','volume3D'}; 
out.File.Load.Loader_auto={'File','Load','Loader_auto'}; 
out.File.Load.Loader_csv={'File','Load','Loader_csv'}; 
out.File.Load.Loader_ioPlugin={'File','Load','Loader_ioPlugin'}; 
out.File.Load.Loader_settings={'File','Load','Loader_settings'}; 
out.File.Load.Loader_sml={'File','Load','Loader_sml'}; 
out.File.Load.Loader_tif={'File','Load','Loader_tif'}; 
out.File.Load.Loader_workflow={'File','Load','Loader_workflow'}; 
out.File.Save.SMLMsaver={'File','Save','SMLMsaver'}; 
out.File.Save.TifSaver={'File','Save','TifSaver'}; 
out.File.Save.saver_csv={'File','Save','saver_csv'}; 
out.File.Save.saver_settings={'File','Save','saver_settings'}; 
out.Process.Assign2C.Get2CIntImagesWF={'Process','Assign2C','Get2CIntImagesWF'}; 
out.Process.Assign2C.Get2CIntLoc={'Process','Assign2C','Get2CIntLoc'}; 
out.Process.Assign2C.Intensity2Channel={'Process','Assign2C','Intensity2Channel'}; 
out.Process.Assign2C.SALMInt2z={'Process','Assign2C','SALMInt2z'}; 
out.Process.Drift.driftcorrection={'Process','Drift','driftcorrection'}; 
out.Process.Drift.showdriftinfo={'Process','Drift','showdriftinfo'}; 
out.Process.Images.Tiff_remove_BG={'Process','Images','Tiff_remove_BG'}; 
out.Process.Modify.Connect2Unconnect={'Process','Modify','Connect2Unconnect'}; 
out.Process.Modify.MathParser={'Process','Modify','MathParser'}; 
out.Process.Modify.Mathematics={'Process','Modify','Mathematics'}; 
out.Process.Modify.Mirror={'Process','Modify','Mirror'}; 
out.Process.Modify.RemoveLocs={'Process','Modify','RemoveLocs'}; 
out.Process.Modify.RenameFields={'Process','Modify','RenameFields'}; 
out.Process.Register.ApplyTransform={'Process','Register','ApplyTransform'}; 
out.Process.Register.RegisterImages={'Process','Register','RegisterImages'}; 
out.Process.Register.RegisterLocs2={'Process','Register','RegisterLocs2'}; 
out.ROIManager.Analyze.AnalyzeCME2Cside={'ROIManager','Analyze','AnalyzeCME2Cside'}; 
out.ROIManager.Analyze.AnalyzeRingCME={'ROIManager','Analyze','AnalyzeRingCME'}; 
out.ROIManager.Analyze.PlotRingCME={'ROIManager','Analyze','PlotRingCME'}; 
out.ROIManager.Analyze.RemoveEmptyCells={'ROIManager','Analyze','RemoveEmptyCells'}; 
out.ROIManager.Analyze.export_tiffs={'ROIManager','Analyze','export_tiffs'}; 
out.ROIManager.Evaluate.CME2CSide={'ROIManager','Evaluate','CME2CSide'}; 
out.ROIManager.Evaluate.CME2DRing={'ROIManager','Evaluate','CME2DRing'}; 
out.ROIManager.Evaluate.CME2DRingFit={'ROIManager','Evaluate','CME2DRingFit'}; 
out.ROIManager.Evaluate.CME3DDSpherefit={'ROIManager','Evaluate','CME3DDSpherefit'}; 
out.ROIManager.Evaluate.CME3DDSpherefit_DC={'ROIManager','Evaluate','CME3DDSpherefit_DC'}; 
out.ROIManager.Evaluate.countingStatistics={'ROIManager','Evaluate','countingStatistics'}; 
out.ROIManager.Evaluate.generalStatistics={'ROIManager','Evaluate','generalStatistics'}; 
out.ROIManager.Segment.segmentCellsCME={'ROIManager','Segment','segmentCellsCME'}; 
out.ROIManager.Segment.segmentNPC={'ROIManager','Segment','segmentNPC'}; 
out.ROIManager.Segment.segmentSitesPolarCME={'ROIManager','Segment','segmentSitesPolarCME'}; 
out.WorkflowModules.Filters.BG_wavelet={'WorkflowModules','Filters','BG_wavelet'}; 
out.WorkflowModules.Filters.DisplayChooser={'WorkflowModules','Filters','DisplayChooser'}; 
out.WorkflowModules.Filters.ImageFilter={'WorkflowModules','Filters','ImageFilter'}; 
out.WorkflowModules.Filters.ImageNormalize={'WorkflowModules','Filters','ImageNormalize'}; 
out.WorkflowModules.Filters.MedianBGcalculator={'WorkflowModules','Filters','MedianBGcalculator'}; 
out.WorkflowModules.Filters.fitOnBackground={'WorkflowModules','Filters','fitOnBackground'}; 
out.WorkflowModules.Fitters.EMCCD_SE_MLE_GPU={'WorkflowModules','Fitters','EMCCD_SE_MLE_GPU'}; 
out.WorkflowModules.Fitters.FalconHD2D={'WorkflowModules','Fitters','FalconHD2D'}; 
out.WorkflowModules.Fitters.RadialSymmetry2D={'WorkflowModules','Fitters','RadialSymmetry2D'}; 
out.WorkflowModules.Fitters.RadialSymmetry3D={'WorkflowModules','Fitters','RadialSymmetry3D'}; 
out.WorkflowModules.Fitters.fitterGUI={'WorkflowModules','Fitters','fitterGUI'}; 
out.WorkflowModules.IntensityCalculator.EvaluateIntensity={'WorkflowModules','IntensityCalculator','EvaluateIntensity'}; 
out.WorkflowModules.IntensityCalculator.IntLoc2pos={'WorkflowModules','IntensityCalculator','IntLoc2pos'}; 
out.WorkflowModules.IntensityCalculator.pushLocs={'WorkflowModules','IntensityCalculator','pushLocs'}; 
out.WorkflowModules.IntensityCalculator.roi2int_fitG={'WorkflowModules','IntensityCalculator','roi2int_fitG'}; 
out.WorkflowModules.IntensityCalculator.roi2int_sumG={'WorkflowModules','IntensityCalculator','roi2int_sumG'}; 
out.WorkflowModules.Loaders.CameraConverter={'WorkflowModules','Loaders','CameraConverter'}; 
out.WorkflowModules.Loaders.GrabFijiStacks={'WorkflowModules','Loaders','GrabFijiStacks'}; 
out.WorkflowModules.Loaders.TifLoader={'WorkflowModules','Loaders','TifLoader'}; 
out.WorkflowModules.Loaders.workflowstarter={'WorkflowModules','Loaders','workflowstarter'}; 
out.WorkflowModules.LocProcessors.LocFilter={'WorkflowModules','LocProcessors','LocFilter'}; 
out.WorkflowModules.LocProcessors.LocSaver={'WorkflowModules','LocProcessors','LocSaver'}; 
out.WorkflowModules.LocProcessors.OnlineReconstruction={'WorkflowModules','LocProcessors','OnlineReconstruction'}; 
out.WorkflowModules.LocProcessors.PlotLocsPreview={'WorkflowModules','LocProcessors','PlotLocsPreview'}; 
out.WorkflowModules.Peakfinders.PeakFinder={'WorkflowModules','Peakfinders','PeakFinder'}; 
out.WorkflowModules.Peakfinders.PeakFinderNMSWF={'WorkflowModules','Peakfinders','PeakFinderNMSWF'}; 
out.WorkflowModules.Peakfinders.RoiAdder={'WorkflowModules','Peakfinders','RoiAdder'}; 
out.WorkflowModules.Peakfinders.RoiCutterWF={'WorkflowModules','Peakfinders','RoiCutterWF'}; 
end 
