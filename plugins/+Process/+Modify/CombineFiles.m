classdef CombineFiles<interfaces.DialogProcessor
    methods
        function obj=CombineFiles(varargin)      
            obj@interfaces.DialogProcessor(varargin{:}) ;  
        end
        
        function out=run(obj,p)
            out=[];
            obj.setPar('undoModule','CombineFiles');
            notify(obj.P,'backup4undo');
            %adjust channel
            fn1=p.file1.Value;
            fn2=p.file2.Value;
            ftarget=min(fn2,fn1);
            fnt=max(fn2,fn1);
            
            ind1=obj.locData.loc.filenumber==fn1;
            ind2=obj.locData.loc.filenumber==fn2;
            if p.addchannel1
                obj.locData.loc.channel(ind1)=obj.locData.loc.channel(ind1)+p.channel1;
            else
                obj.locData.loc.channel(ind1)=p.channel1;
            end
            if p.addchannel2
                obj.locData.loc.channel(ind2)=obj.locData.loc.channel(ind2)+p.channel2;
            else
                obj.locData.loc.channel(ind2)=p.channel2;
            end
            
            % adjust filenumber
            obj.locData.loc.filenumber(ind2)=ftarget;
            %rename file
            obj.locData.files.file(ftarget).name=strrep(obj.locData.files.file(ftarget).name,p.file1.selection,p.newfilename);
            %adjust file info (e.g. number of frames, ROI etc), combine tif
            obj.locData.files.file(ftarget).tif=[obj.locData.files.file(ftarget).tif obj.locData.files.file(fnt).tif];
             obj.locData.files.file(ftarget).numberOfTif=length(obj.locData.files.file(ftarget).tif);
            obj.locData.files.file(ftarget).raw=[obj.locData.files.file(ftarget).raw obj.locData.files.file(fnt).raw];
            %etc, if different pixelsize: warn
            %remove second file.
            obj.locData.files.file(fnt)=[];
            obj.locData.files.filenumberEnd=obj.locData.files.filenumberEnd-1;
            %update filelist
            obj.locData.regroup;
            initGuiAfterLoad(obj)
             
        end
        function pard=guidef(obj)
            pard=guidef(obj);
        end
    end
end

function newfilenameb_callback(a,b,obj)
p=obj.getGuiParameters;
f1=p.file1.selection;
f2=p.file2.selection;
re='[^_]*';
[i1,i2]=regexp(f1,re);
for k=1:length(i1)
    f2=regexprep(f2,f1(i1(k):i2(k)),'');
end

fo=['C' num2str(p.channel2) '_' f2 'C' num2str(p.channel1) '_' f1];
fo=regexprep(fo,'[_]*','_');
po.newfilename=fo;
obj.setGuiParameters(po)

end

function pard=guidef(obj)
pard.t1.object=struct('String','Channel','Style','text');
pard.t1.position=[1,4];
pard.t1.Width=0.5;

pard.t3.object=struct('String','add ','Style','text');
pard.t3.position=[1,4.5];
pard.t3.Width=0.5;

pard.file1.object=struct('Style','popupmenu','String','File');
pard.file1.position=[2,1];
pard.file1.object.TooltipString='choose localization file data set';
pard.file1.Width=3;

pard.file2.object=struct('Style','popupmenu','String','File');
pard.file2.position=[3,1];
pard.file2.object.TooltipString='choose localization file data set';
pard.file2.Width=3;


pard.channel1.object=struct('String','1','Style','edit');
pard.channel1.position=[2,4];
pard.channel1.Width=0.5;

pard.channel2.object=struct('String','2','Style','edit');
pard.channel2.position=[3,4];
pard.channel2.Width=0.5;

pard.addchannel1.object=struct('String','','Style','checkbox');
pard.addchannel1.position=[2,4.5];
pard.addchannel1.Width=0.5;

pard.addchannel2.object=struct('String','','Style','checkbox');
pard.addchannel2.position=[3,4.5];
pard.addchannel2.Width=0.5;

pard.newfilenameb.object=struct('String','New filename:','Style','pushbutton','Callback',{{@newfilenameb_callback,obj}});
pard.newfilenameb.position=[5,1];

pard.newfilename.object=struct('String','newfile_sml','Style','edit');
pard.newfilename.position=[5,2];
pard.newfilename.Width=3;

pard.syncParameters={{'filelist_short','file1',{'String'}},{'filelist_short','file2',{'String'}}};

pard.plugininfo.type='ProcessorPlugin';
end