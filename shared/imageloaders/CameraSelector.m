classdef CameraSelector<handle
    %UNTITLED2 Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        handle
        imloader
    end
    
    methods
        function obj=CameraSelector
            makeGui(obj)
        end
    end
    
end

function makeGui(obj)
width=600;
height=500;
lineheight=25;
posbutton=width*.8;
buttonwidth=width*.15;
if isempty(obj.handle)||~isvalid(obj.handle)
     obj.handle=figure('Units','normalized','Units','pixels','Position',[150,200,width,height],'Name','CameraSelector','NumberTitle','off');
     delete(obj.handle.Children);
end
%load images
hp=uicontrol('Style','pushbutton','String','Load images','Position',[posbutton height-100,buttonwidth,lineheight],'Callback',{@loadimages,obj});

tcam=uitable(obj.handle,'Position',[10 height-lineheight*3 width-200,lineheight*2]);
tcam.ColumnName={'Camera Name','ID field','ID'};
tcam.Data={'Cam1','Cam_ID','001'};
tcam.CellSelectionCallback={@cellselect,obj};
end
function cellselect(table,data,obj)
if data.Indices(2)==2
    md=obj.imloader.metadata.allmetadata;
    fo=browsefields(md);
end

end

function loadimages(a,b,obj)
[file path]=uigetfile('*.*');
if file
obj.imloader=imageloaderAll([path file]);
end
end