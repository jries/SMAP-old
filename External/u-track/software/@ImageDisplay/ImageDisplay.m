classdef ImageDisplay < MovieDataDisplay
    %Abstract class for displaying image processing output
%
% Copyright (C) 2016, Danuser Lab - UTSouthwestern 
%
% This file is part of u-track.
% 
% u-track is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
% 
% u-track is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
% 
% You should have received a copy of the GNU General Public License
% along with u-track.  If not, see <http://www.gnu.org/licenses/>.
% 
% 
    properties
        Colormap ='gray';
        Colorbar ='off';
        ColorbarLocation ='EastOutside';
        CLim = [];
        Units='';
        sfont = {'FontName', 'Helvetica', 'FontSize', 18};
        lfont = {'FontName', 'Helvetica', 'FontSize', 22};
        ScaleFactor = 1;
        NaNColor = [0 0 0];
    end
    methods
        function obj=ImageDisplay(varargin)
            obj@MovieDataDisplay(varargin{:});
        end
        
        function h=initDraw(obj,data,tag,varargin)
            % Plot the image and associate the tag
            h=imshow(data/obj.ScaleFactor,varargin{:});
            set(h,'Tag',tag,'CDataMapping','scaled');
            hAxes = get(h,'Parent');
            set(hAxes,'XLim',[0 size(data,2)],'YLim',[0 size(data,1)]);
            
            % Tag all objects displayed by this class, so they can be
            % easily identified and cleared.
            userData = get(h,'UserData');
            userData.DisplayClass = 'ImageDisplay';
            set(h,'UserData',userData)
            
            % Clean existing image and set image at the bottom of the stack
            child = get(hAxes,'Children');
            imChild = child(strcmp(get(child,'Type'),'image'));
            delete(imChild(imChild ~= h));
            uistack(h, 'bottom');
            
            obj.applyImageOptions(h)
        end
        
        function updateDraw(obj,h,data)
            if(obj.ScaleFactor ~= 1)
                set(h,'CData',data/obj.ScaleFactor)
            else
                set(h,'CData',data);
            end
            obj.applyImageOptions(h)
        end
        
        function applyImageOptions(obj,h)
            hAxes = get(h,'Parent');
            % Set the colormap
            if any(isnan(get(h, 'CData')))
                c = colormap(obj.Colormap);
                c=[obj.NaNColor; c];
                colormap(hAxes, c);
            else
                colormap(hAxes,obj.Colormap);
            end
            
            % Set the colorbar
            hCbar = findobj(get(hAxes,'Parent'),'Tag','Colorbar');
            axesPosition = [0 0 1 1];
            if strcmp(obj.Colorbar,'on')
                if length(obj.ColorbarLocation) >6 && strcmp(obj.ColorbarLocation(end-6:end),'Outside'),
                    axesPosition = [0.05 0.05 .9 .9];
                end
                if isempty(hCbar)
                    set(hAxes,'Position',axesPosition);
                    hCbar = colorbar('peer',hAxes,obj.sfont{:});
                end
                set(hCbar,'Location',obj.ColorbarLocation);
                ylabel(hCbar,obj.Units,obj.lfont{:});
            else
                if ~isempty(hCbar),colorbar(hCbar,'delete'); end
                set(hAxes,'Position',axesPosition);
            end
            
            % Set the color limits
            if ~isempty(obj.CLim),
                set(hAxes,'CLim',obj.CLim/obj.ScaleFactor);
            else
                set(hAxes,'CLim', [0 1]);
            end
        end
    end
    
    methods (Static)
        function params=getParamValidators()
            params(1).name='Colormap';
            params(1).validator=@ischar;
            params(2).name='Colorbar';
            params(2).validator=@(x) any(strcmp(x,{'on','off'}));
            params(3).name='CLim';
            params(3).validator=@isvector;
            params(4).name='Units';
            params(4).validator=@ischar;
            params(5).name='sfont';
            params(5).validator=@iscell;
            params(6).name='lfont';
            params(6).validator=@iscell;
            params(7).name='ColorbarLocation';
            params(7).validator=@(x) any(strcmp(x, ImageDisplay.getColorBarLocations()));
            params(8).name='ScaleFactor';
            params(8).validator=@isscalar;
            params(9).name='NaNColor';
            params(9).validator=@isvector;
        end
        
        function locations = getColorBarLocations()
            locations = {...
                'North', 'South', 'East', 'West',...
                'NorthOutside', 'SouthOutside', 'EastOutside',...
                'WestOutside' };
        end
        
        function f=getDataValidator()
            f=@isnumeric;
        end
    end
end