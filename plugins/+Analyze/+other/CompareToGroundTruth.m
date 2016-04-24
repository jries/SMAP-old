classdef CompareToGroundTruth<interfaces.DialogProcessor
    methods
        function obj=CompareToGroundTruth(varargin)        
                obj@interfaces.DialogProcessor(varargin{:});
            obj.inputParameters={'numberOfLayers','sr_layerson','mainfile','cam_pixelsize_nm'};
%             end   
        end
        
        function out=run(obj,p)
            out=[];
            %get localizations
            [path,file,ext]=fileparts(p.mainfile);
            if ~exist(path,'dir')
                path='settings';
            end
            filenew=fullfile(path,[file '_temp.csv']);
            
            lochere=obj.locData.copy;
            if p.shiftpix
            shiftx=-0.5*p.cam_pixelsize_nm;
            shifty=-0.5*p.cam_pixelsize_nm;
            else
                shiftx=0;
                shifty=0;
            end
            lochere.loc.xnm=lochere.loc.xnm+shiftx;
            lochere.loc.ynm=lochere.loc.ynm+shifty;
            descfile=saveLocalizationsCSV(lochere,filenew,p.onlyfiltered,p.numberOfLayers,p.sr_layerson);
            
            %modify challenge data
            challengepath=['External' filesep 'SMLMChallenge' filesep];
            javapath=['"' pwd filesep challengepath 'challenge.jar"'];

            settingsfile=[challengepath 'CompareLocalizationSettings.txt'];
            
            replacements={'firstRow1','0','shiftY1','0','txtFile1',strrep(filenew,'\','/'),'colY1','2','shiftX1','0','colF1','0','colI1','4','shiftUnit1','nm','txtDesc1',strrep(descfile,'\','/')};
            modifySettingsFile(settingsfile,replacements{:});
            oldp=pwd;
            cd(challengepath);
            disp('to contintue with Matlab, close SMLMChallenge application');
            system(['java -jar ' javapath]) 
            %later fix jave program and call via
            %smlm.assessment.application.Application
            %after adding javaclasspath(javapath)
            cd(oldp)
        end
        function pard=guidef(obj)
            pard=guidef;
        end
    end
end




function pard=guidef
pard.onlyfiltered.object=struct('Style','checkbox','String','Export filtered (displayed) localizations.','Value',1);
pard.onlyfiltered.position=[2,1];
pard.onlyfiltered.Width=4;

pard.shiftpix.object=struct('Style','checkbox','String','Shift by 0.5 camera pixels','Value',0);
pard.shiftpix.position=[3,1];
pard.shiftpix.Width=4;
end