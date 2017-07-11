classdef calibrate3D_GUI<handle
    properties
        guihandles
    end
    methods
        function obj=calibrate3D_GUI(varargin)  
            %constructur: make GUI
            addpath('shared')
            h=figure('Name','3D calibration','MenuBar','none','ToolBar','none');
            h.Position(3:4)=[450, 600];
            top=h.Position(4);
            vsep=24;
            if ispc
                fontsize=12;
            else 
                fontsize=14;
            end
            xpos1=10;
            xw=100;
            obj.guihandles.title=uicontrol('style','text','String','Calibrate PSF model for MLE fit from bead stacks. ','Position',[xpos1,top-vsep,xw*3,vsep],'FontSize',fontsize);
            
            obj.guihandles.selectfiles=uicontrol('style','pushbutton','String','Select camera files','Position',[xpos1,top-2*vsep,xw*1.5,vsep],'FontSize',fontsize,'Callback',@obj.selectfiles_callback);
            obj.guihandles.filelist=uicontrol('style','listbox','String','','Position',[xpos1+1.5*xw,top-4*vsep,xw*2.5,vsep*3],'FontSize',fontsize);
            obj.guihandles.selectoutputfile=uicontrol('style','pushbutton','String','Select otuput file','Position',[xpos1,top-5*vsep,xw*1.5,vsep],'FontSize',fontsize,'Callback',@obj.selectoutputfile_callback);
            obj.guihandles.outputfile=uicontrol('style','edit','String','data/bead_3dcal.mat','Position',[xpos1+1.5*xw,top-5*vsep,xw*2.5,vsep],'FontSize',fontsize);
            
            hatitle='left';
            ha='right';
            obj.guihandles.csplinet=uicontrol('style','text','String','General paramters: ','Position',[xpos1,top-7*vsep,xw*4,vsep],'FontSize',fontsize,'HorizontalAlignment',hatitle,'FontWeight','bold');
            obj.guihandles.dzt=uicontrol('style','text','String','Distance between frames (nm)','Position',[xpos1,top-8*vsep,xw*2.5,vsep],'FontSize',fontsize,'HorizontalAlignment',ha);
            obj.guihandles.dz=uicontrol('style','edit','String','10','Position',[xpos1+2.5*xw,top-8*vsep,xw*1,vsep],'FontSize',fontsize);
            
            obj.guihandles.modalityt=uicontrol('style','text','String','3D modality ','Position',[xpos1,top-9*vsep,xw*2.5,vsep],'FontSize',fontsize,'HorizontalAlignment',ha);
            obj.guihandles.modality=uicontrol('style','popupmenu','String',{'astigmatism','arbitrary','2D PSF'},'Position',[xpos1+2.5*xw,top-9*vsep,xw*1.5,vsep],'FontSize',fontsize,'Callback',@obj.modality_callback);
            
            obj.guihandles.corrzt=uicontrol('style','text','String','Correct bead z-positions using ','Position',[xpos1,top-10*vsep,xw*2.5,vsep],'FontSize',fontsize,'HorizontalAlignment',ha);
            obj.guihandles.corrzselect=uicontrol('style','popupmenu','String',{'none','cross-correlation','shape (astig)'},...
                'Value',3,'Position',[xpos1+2.5*xw,top-10*vsep,xw*1.5,vsep],'FontSize',fontsize,'Callback',@obj.zcorr_callback);
            
            obj.guihandles.zcorrframest=uicontrol('style','text','String','frames to use for CC: ','Position',[xpos1+1.5*xw,top-11*vsep,xw*2,vsep],'FontSize',fontsize,'Visible','off','HorizontalAlignment',ha);
            obj.guihandles.zcorrframes=uicontrol('style','edit','String','30','Position',[xpos1+3.5*xw,top-11*vsep,xw*.5,vsep],'FontSize',fontsize,'Visible','off');
            
            obj.guihandles.filtert=uicontrol('style','text','String','Filter size for peak finding','Position',[xpos1,top-12*vsep,xw*2.5,vsep],'FontSize',fontsize,'HorizontalAlignment',ha);
            obj.guihandles.filter=uicontrol('style','edit','String','2','Position',[xpos1+2.5*xw,top-12*vsep,xw*1,vsep],'FontSize',fontsize);
            
            
     
            obj.guihandles.csplinet=uicontrol('style','text','String','Cspline paramters: ','Position',[xpos1,top-14*vsep,xw*4,vsep],'FontSize',fontsize,'HorizontalAlignment',hatitle,'FontWeight','bold');
            obj.guihandles.roisizet=uicontrol('style','text','String','ROI size: X,Y (pixels): ','Position',[xpos1,top-15*vsep,xw*2,vsep],'FontSize',fontsize,'HorizontalAlignment',ha);
            obj.guihandles.ROIxy=uicontrol('style','edit','String','21','Position',[xpos1+2*xw,top-15*vsep,xw*.5,vsep],'FontSize',fontsize);
            obj.guihandles.roisizezt=uicontrol('style','text','String','Z (frames): ','Position',[xpos1+2.5*xw,top-15*vsep,xw,vsep],'FontSize',fontsize,'HorizontalAlignment',ha);
            obj.guihandles.ROIz=uicontrol('style','edit','String','201','Position',[xpos1+3.5*xw,top-15*vsep,xw*.5,vsep],'FontSize',fontsize);
            
            obj.guihandles.smootht=uicontrol('style','text','String','Smoothing parameter in Z: ','Position',[xpos1,top-16*vsep,xw*2,vsep],'FontSize',fontsize,'HorizontalAlignment',ha);
%             obj.guihandles.smoothxy=uicontrol('style','edit','String','0','Position',[xpos1+2*xw,top-16*vsep,xw*.5,vsep],'FontSize',fontsize);
%             obj.guihandles.smoothzt=uicontrol('style','text','String','Z: ','Position',[xpos1+2.5*xw,top-16*vsep,xw,vsep],'FontSize',fontsize,'HorizontalAlignment',ha);
%             obj.guihandles.smoothz=uicontrol('style','edit','String','2','Position',[xpos1+3.5*xw,top-16*vsep,xw*.5,vsep],'FontSize',fontsize);
            obj.guihandles.smoothz=uicontrol('style','edit','String','2','Position',[xpos1+2*xw,top-16*vsep,xw*.5,vsep],'FontSize',fontsize);
            
            obj.guihandles.gausst=uicontrol('style','text','String','Gauss fit parameters: ','Position',[xpos1,top-18*vsep,xw*4,vsep],'FontSize',fontsize,'HorizontalAlignment',hatitle,'FontWeight','bold');
            obj.guihandles.gaussmint=uicontrol('style','text','String','Range (nm). minimum: ','Position',[xpos1,top-19*vsep,xw*2,vsep],'FontSize',fontsize,'HorizontalAlignment',ha);
            obj.guihandles.gaussmin=uicontrol('style','edit','String','-500','Position',[xpos1+2*xw,top-19*vsep,xw*.5,vsep],'FontSize',fontsize);
            obj.guihandles.gaussmaxt=uicontrol('style','text','String','maximum: ','Position',[xpos1+2.5*xw,top-19*vsep,xw,vsep],'FontSize',fontsize,'HorizontalAlignment',ha);
            obj.guihandles.gaussmax=uicontrol('style','edit','String','500','Position',[xpos1+3.5*xw,top-19*vsep,xw*.5,vsep],'FontSize',fontsize);
             obj.guihandles.gaussroit=uicontrol('style','text','String','ROI size (pixels): ','Position',[xpos1,top-20*vsep,xw*2,vsep],'FontSize',fontsize,'HorizontalAlignment',ha);
            obj.guihandles.gaussroi=uicontrol('style','edit','String','17','Position',[xpos1+2*xw,top-20*vsep,xw*.5,vsep],'FontSize',fontsize);
           
            obj.guihandles.run=uicontrol('style','pushbutton','String','Calculate bead calibration','Position',[xpos1,top-22*vsep,xw*4,vsep],'FontSize',fontsize,'Callback',@obj.run_callback);
            obj.guihandles.help=uicontrol('style','pushbutton','String','Help','Position',[xpos1+xw,top-23*vsep,xw*2,vsep],'FontSize',fontsize,'Callback',@obj.help_callback);
                      
            obj.guihandles.status=uicontrol('style','text','String','Status','Position',[xpos1,top-25*vsep,xw*4,vsep],'FontSize',fontsize,'HorizontalAlignment','left');
            
        end
        function selectfiles_callback(obj,a,b)
            sf=selectManyFiles;
            sf.guihandles.filelist.String=(obj.guihandles.filelist.String);
            waitfor(sf.handle);
            obj.guihandles.filelist.String=sf.filelist;
            obj.guihandles.filelist.Value=1;
            if isempty(obj.guihandles.outputfile.String)
                [path,file]=fileparts(sf.filelist{1});
                obj.guihandles.outputfile.String=[path file '_3Dcorr.mat'];
                
            end
        end
        function selectoutputfile_callback(obj,a,b)
            of=obj.guihandles.outputfile.String;
            if isempty(of)
                of='_corr3D.mat';
            end
            [f,p]=uiputfile(of);
            if f
            obj.guihandles.outputfile.String=[p,f];
            end
        end
        function modality_callback(obj,a,b)
            astigonly={'gausst','gaussmin','gaussmint','gaussmax','gaussmaxt','gaussroi','gaussroit'};
            switch obj.guihandles.modality.String{obj.guihandles.modality.Value}
                case 'astigmatism'
                    vis='on';
                    obj.guihandles.corrzselect.String={'none','cross-correlation','shape (astig)'};
                    
                otherwise
                    vis='off';
                    obj.guihandles.corrzselect.String={'none','cross-correlation'};
                    obj.guihandles.corrzselect.Value=min(obj.guihandles.corrzselect.Value,2);
            end
            for k=1:length(astigonly)
                obj.guihandles.(astigonly{k}).Visible=vis;
            end
        end
        function zcorr_callback(obj,a,b)
            corrz={'zcorrframest','zcorrframes'};
            if contains(obj.guihandles.corrzselect.String{obj.guihandles.corrzselect.Value},'cross')
                vis='on';
            else
                vis='off';
            end
            for k=1:length(corrz)
                obj.guihandles.(corrz{k}).Visible=vis;
            end
        end
        function out=run_callback(obj,a,b)
            p.filelist=obj.guihandles.filelist.String;
            p.outputfile=obj.guihandles.outputfile.String;
            p.dz=str2double(obj.guihandles.dz.String);
            p.modality=obj.guihandles.modality.String{obj.guihandles.modality.Value};
            p.zcorr=obj.guihandles.corrzselect.String{obj.guihandles.corrzselect.Value};
            p.ROIxy=str2double(obj.guihandles.ROIxy.String);
            p.ROIz=str2double(obj.guihandles.ROIz.String);
%             p.smoothxy=str2double(obj.guihandles.smoothxy.String);
            p.smoothxy=0;
            p.smoothz=str2double(obj.guihandles.smoothz.String);
            p.gaussrange=[str2double(obj.guihandles.gaussmin.String) str2double(obj.guihandles.gaussmax.String)];
            p.filtersize=str2double(obj.guihandles.filter.String);
            p.zcorrframes=str2double(obj.guihandles.zcorrframes.String);
            p.gaussroi=str2double(obj.guihandles.gaussroi.String);
            p.status=obj.guihandles.status;
            if isempty(p.filelist)
                warndlg('please select image files first')
                return
            end
            
            calibrate3D(p);
        end    
        function help_callback(obj,a,b)
            helpstring={'write documentation'};
            f=figure;
            h=uicontrol('Parent',f,'Style','text','Units','normalized','Position',[0 0 1 1],'HorizontalAlignment','Left');
            h.String=textwrap(h,helpstring);
            
        end
    end
end


