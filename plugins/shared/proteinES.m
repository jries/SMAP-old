classdef proteinES < handle
    % pes = proteinES('121212', 'sla1', bu, 'ellipse','bar')
    % pes.addTimePoint(table([1],[1],[10],[20],[5],[5],[200],{'bar'},'VariableNames',{'time', 'dz', 'th', 'de', 'zm', 'shi', 'nm', 'sha'}));
    properties
        %define class properties if needed
        proteinID
        proteinName
        endocyticPM
        timePar
        capShape
        bottomShape
        timeData
    end
    methods
        function obj=proteinES(proteinID, proteinName, basicUnit, capShape, bottomShape)   %replace by filename        
            obj.proteinID = proteinID;
            obj.proteinName = proteinName;
            obj.endocyticPM = basicUnit;
            obj.capShape = capShape;
            obj.bottomShape = bottomShape;
            timePar = table([],[],[],[],[],[],[],[],'VariableNames',{'time', 'dz', 'th', 'de', 'zm', 'shi', 'nm', 'sha'});
            % obj.timePar = table([1],[1],[10],[20],[5],[5],[200],{'bar'},'VariableNames',{'time', 'dz', 'th', 'de', 'zm', 'shi', 'nm', 'sha'});
            obj.timePar = timePar;
            obj.timeData = [];
        end
        function addTimePoint(obj, aTable)
            obj.timePar = [obj.timePar;aTable];
        end
        function img = getTimeImg(obj, t, viewType, imageSize)
            if isempty(obj.timeData.(['t' int2str(t)])) 
                obj.timeData.(['t' t]) = obj.getTime(t);
            end    
            switch viewType
                case 'top'
                    % accumulate all the even points in the space to z=1
                    idx = int16([obj.timeData.(['t' int2str(t)]).x obj.timeData.(['t' int2str(t)]).y])+imageSize/2;
                    idxZ = obj.timeData.(['t' int2str(t)]).z;
                    proImg = accumarray(idx, idxZ,[],@(x) length(x));
                    proImg = [proImg zeros(min(size(proImg)),imageSize-min(size(proImg)))];
                    img = [proImg; zeros(imageSize-min(size(proImg)),imageSize)];

                    figure(31)
                    image(img, 'CDataMapping', 'scaled')
                    % imwrite(uint8(round(mat2gray(proImg)*255)), [p.folderPath '\' p.imgPath])
                case 'side'
                    % accumulate all the even points in the space to z=1
                    idx = int16([obj.timeData.(['t' int2str(t)]).x+imageSize/2 obj.timeData.(['t' int2str(t)]).z+1]);
                    idxZ = obj.timeData.(['t' int2str(t)]).y+imageSize/2;
                    proImg = accumarray(idx, idxZ,[],@(x) length(x));
                    Size = size(proImg);
                    proImg = [proImg; zeros(imageSize-Size(1),Size(2))];
                    Size = size(proImg);
                    img = [proImg zeros(imageSize, imageSize-Size(2))];
                    figure(31)
                    image(img', 'CDataMapping', 'scaled')
                    % imwrite(uint8(round(mat2gray(proImg')*255)), [p.folderPath '\' p.imgPath])
            end
        end
        function time = getTime(obj, t)
            % timeTime = getTime(pes,1);
            thisTP = obj.timePar(obj.timePar.time==t,:);
            inner = basicUnit(obj.endocyticPM.capDepth+thisTP.shi, obj.endocyticPM.capShape, obj.endocyticPM.bottomDepth, obj.endocyticPM.bottomShape, obj.endocyticPM.radius+thisTP.shi);
            
            if strcmp(obj.endocyticPM.capShape, obj.capShape) && strcmp(obj.endocyticPM.bottomShape, obj.bottomShape) % Check the shapes of inner and outer columns are the same or not
                outer = basicUnit(obj.endocyticPM.capDepth+thisTP.shi+thisTP.th, obj.endocyticPM.capShape, obj.endocyticPM.bottomDepth, obj.endocyticPM.bottomShape, obj.endocyticPM.radius+thisTP.shi+thisTP.th);
            end
            
            Px = fullPixel(outer);
            inOuter = basicUnit.pointsInCy(outer.basicComponent3D.x, outer.basicComponent3D.y, outer.basicComponent3D.z, Px.x(:), Px.y(:), Px.z(:));
            Px = Px(inOuter,:);
            inInner = basicUnit.pointsInCy(inner.basicComponent3D.x, inner.basicComponent3D.y, inner.basicComponent3D.z, Px.x(:), Px.y(:), Px.z(:));
            Px = Px(~inInner,:);
            thisBasicComponent3DPx = Px;
            time = thisBasicComponent3DPx(thisBasicComponent3DPx.z <= thisTP.zm+thisTP.de/2 & thisBasicComponent3DPx.z >= thisTP.zm-thisTP.de/2 ,:);
        end
    end
end
