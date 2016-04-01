classdef Undo< interfaces.GuiModuleInterface & interfaces.LocDataInterface
    properties
        locDataOld
        undoModule
    end
    methods
        function obj=Undo(varargin)
            obj@interfaces.GuiModuleInterface(varargin{:})          
        end
        function makeGui(obj)
            obj.guihandles.undobutton=uicontrol(obj.handle,'Style','pushbutton','String','Undo','Units','normalized',...
                'Position',[0.8,0.005,.07,.03],'Callback',@obj.undo_callback);
             obj.guihandles.undobutton.Units='pixels';
             obj.guihandles.undobutton.Position(4)=28;
            addlistener(obj.P,'backup4undo',@obj.backup);
        end
        function backup(obj,event,data) 
            obj.undoModule=obj.getPar('undoModule');
            obj.locDataOld=obj.locData.copy;
        end
        function undo_callback(obj,a,b)
            if isempty(obj.locDataOld)
%                 disp('nothing stored for undo')
                obj.status(['nothing stored for undo ' ])
                return
            end        
            temp=obj.locData.copy;
            obj.locData.loc=obj.locDataOld.loc;
            obj.locData.grouploc=obj.locDataOld.grouploc;
            obj.locData.files=obj.locDataOld.files;
            obj.setPar('locFields',fieldnames(obj.locData.loc));
            obj.locDataOld=temp;
            
            obj.status(['undo performed: ' obj.undoModule])
            obj.undoModule=['undo of ' obj.undoModule];
            
        end
    end
end

