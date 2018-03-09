classdef boundaryFinder<interfaces.SEEvaluationProcessor
    properties
        boundary
    end
    methods
        function obj=boundaryFinder(varargin)        
                obj@interfaces.SEEvaluationProcessor(varargin{:});
        end
        function out=run(obj, inp)
            % Import the kimograph
            test = obj.site.evaluation.BALM_fibril_growth.kimograph;
            test(test>=1)=1;

            slideStep = inp.slideStep;
            gridRange = inp.gridRange;
            adjM = inp.adjM;

            % expand test
            oriSize = size(test);
            expand = zeros(oriSize(1),gridRange*slideStep);
            test = [test, expand];

            % Binarize the kimograph
            [testRefx, testRefy] = find(test);
            [testNZx, testNZy] = find(test);

            % use euclidean distances to filter out noises
            K = round(length(testRefx)*5/100);

            ED = pdist2([testNZx testNZy],[testNZx testNZy],'euclidean');
            sortedED = sort(ED,2);
            sortedED = sortedED(:,1:K);
            rmean = mean(sortedED,2);

            EDcutoff = quantile(rmean, 0.95);

            test(sub2ind(size(test), testNZx([rmean>=EDcutoff]), testNZy([rmean>=EDcutoff]))) = 0;

            % Get cordinates of the kimograph
            Size = size(test);
            [CorX CorY] = meshgrid(1:Size(1), 1:Size(2));
            Order = sub2ind(Size, CorX(:), CorY(:));

            % Make grids (see gridRange)
            adjCorX = ceil(CorX(:)./gridRange);
            adjCorY = ceil(CorY(:)./gridRange);

            % Get density
            D = accumarray([adjCorX adjCorY], test(Order), [], @sum);

            % Filter out grids that might be background
            cutoff = 3; %#Par
            D(D<=cutoff)=0;

            % Do sliding subtraction (see slideStep)
            % to see the next slideStep-1 right neighbor is higher or not than the
            % target
            SizeD = size(D);

            diffOne = cell(slideStep-1,1);
            diffAll = zeros(SizeD(1), SizeD(2)-slideStep);
            for i = 1:(slideStep-1)
                diffAll = diffAll + (D(:,1+i:end-slideStep+i) - D(:,1:(end-slideStep)) < 0);
            end

            % Get grids with 4 zero grids on its right-hand side
            diffAllLog = diffAll == slideStep-1; 
            [bX,bY] = find(diffAllLog, 1, 'last');

            % Get boundary (itself is 1, and its cumsum is also 1)
            boundary = diffAllLog & cumsum(diffAllLog,2,'reverse') == 1;

            %figure(50)
            %mesh(boundary)
            %rotate3d on

            % Find the corrdinates of boundary
            [Cx, Cy]= find(boundary,sum(boundary(:)),'last');

            % Make sure every frame (x-axies) have a point
            bounSize = size(boundary);
            CxFinal = (1:bounSize(1))';
            CyFinal = zeros(bounSize(1),1);
            CyFinal(Cx) = Cy;

            % Find the cummax (next points should always >= the previous points)
            CyFinal = cummax(CyFinal);
            CxGridParent = CxFinal*gridRange;
            CyGridParent = CyFinal*gridRange;

            meaningfulX = (CxGridParent+1)<=Size(1);
            CxGridParent = CxGridParent(meaningfulX);
            CyGridParent = CyGridParent(meaningfulX);

            CxFParent = (0:Size(1)-1)';
            CyFParent = zeros(Size(1),1);
            CyFParent(CxGridParent+1) = CyGridParent;
            CyFParent = cummax(CyFParent);


            %figure(50)
            %scatter(testNZx, testNZy)
            %hold on
            %plot(CxFParent, CyFParent)


            Mref = calMeasurement(CxFParent,CyFParent,testNZx,testNZy,Size);
            MrefO = Mref;
            CyFParentA = CyFParent;
            for i=size(CxFParent):-1:1
               thisFrameY = testNZy(testNZx == CxFParent(i));
               thisFrameY = thisFrameY(thisFrameY>CyFParent(i));
               ii = 1;
               futherStep = 0;
               while ii <= length(thisFrameY)
                   CyFParentO = CyFParentA;

                   % New boundary
                   CyFParentA(i)=thisFrameY(ii);

                   CyFParentA = cummax(CyFParentA);

                   % fit to new boundary
                   Mnew = calMeasurement(CxFParent,CyFParentA,testNZx,testNZy,Size);
                   if Mnew >= MrefO*adjM && Mnew >= Mref %#Par
                       MrefO = Mnew;
                       futherStep = 0;
                   else
                       if futherStep <= 5
                        CyFParentA = CyFParentO;
                        Mnew = MrefO;
                        futherStep = futherStep+1;
                       else
                        CyFParentA = CyFParentO;
                        Mnew = MrefO;
                        break
                       end
                   end
                   ii=ii+1;
               end
                if isinteger(i/10)
                    figure(51)
                    scatter(testNZx, testNZy)
                    hold on
                    plot(CxFParent, CyFParentA)
                    hold off
                end
            end
            
            KimoPar = obj.site.evaluation.BALM_fibril_growth;
            fig = KimoPar.kimograph;
            co=quantile(fig(:),0.999);
            fig(fig>co)=co;
            out.boundary = [CxFParent CyFParentA];

            % Remove redundant points (in stalled region)
            stepWidth = diff(out.boundary(:,2));
            indStepWidthNon0 = find(stepWidth);
            indPoint2Keep = [indStepWidthNon0; indStepWidthNon0+1];
            indPoint2Keep = unique(indPoint2Keep);
            out.boundary = out.boundary([1; indPoint2Keep; end],:);
            
            % Merge steps by move
            while 1
                [newboundary, idx] =  mergeStallsByMove(out.boundary);
                idxRef = idx;
                if isequal(idxRef, idx)
                    break
                end
                out.boundary =  newboundary;
            end
            
            % Merge steps by move
            while 1
                [newboundary, idx] =  mergeByMove(out.boundary);
                idxRef = idx;
                if isequal(idxRef, idx)
                    break
                end
                out.boundary =  newboundary;
            end
            
            % Merge steps by time
            out.boundary =  mergeByTime(out.boundary);
            
            % Merge steps by move
            out.boundary = mergeStallsByMove(out.boundary);
            
            % Merge steps by move again
            while 1
                [newboundary, idx] =  mergeByMove(out.boundary);
                idxRef = idx;
                if isequal(idxRef, idx)
                    break
                end
                out.boundary =  newboundary;
            end
            
            % Merge steps by time
            out.boundary =  mergeByTime(out.boundary);
            
            % Remove redundant points (in stalled region)
            stepWidth = diff(out.boundary(:,2));
            indStepWidthNon0 = find(stepWidth);
            indPoint2Keep = [indStepWidthNon0; indStepWidthNon0+1];
            indPoint2Keep = unique(indPoint2Keep);
            out.boundary = out.boundary([1; indPoint2Keep; end],:);
            
            out.stallTime = calStallTime(out.boundary);
            out.stepWidth = calStepWidth(out.boundary);
            out.growthTime = calGrowthTime(out.boundary);
            out.avgRate = calAvgRate(out.boundary);
            
            
            h=obj.setoutput('kimograph');
            imagesc(h,(fig))
            hold(h,'on')
            plot(h,CyFParentA,CxFParent, 'LineWidth', 1, 'Color', 'g')
            plot(h,out.boundary(:,2),out.boundary(:,1), 'LineWidth', 1.5, 'Color', 'w')
            hold(h,'off')            
            

            %xlabel(h,'xnm')
            %ylabel(h,'frame')
            %xticklabels(h, KimoPar.xn)
            %yticklabels(h, KimoPar.fr)


            
            if 0
                out.avgRate = (CyFParentA(end)-CyFParentA(1))/(CxFParent(end)-CxFParent(1));
                out.stepWidth = stepWidthNon0;
                out.stallTime = stallTime;

                h2=obj.setoutput('statistics');
                axes(h2);
                ax1 = subplot(1,2,1);
                histogram(ax1, stepWidthNon0, 'BinWidth', 1);
                title(ax1,'Rate');
                ax2 = subplot(1,2,2);
                histogram(ax2, stallTime, 'BinWidth', 1);
                title(ax2,'Stall time');
            end
        end
     
        function pard=guidef(obj)
            pard=guidef;
        end
    end

end



function pard=guidef
pard.t_gridRange.object=struct('Style','text','String','Bin size of the grids');
pard.t_gridRange.position=[1,1];
pard.t_gridRange.Width=2;

pard.gridRange.object=struct('Style','edit','String',5);
pard.gridRange.position=[1,3];
pard.gridRange.TooltipString = 'If you set it as 5, it means before the density comparison, every grid will be set to cover a 5-by-5 area in the original coordinates';
            
pard.t_slideStep.object=struct('Style','text','String','Slide step(s)');
pard.t_slideStep.position=[2,1];
pard.t_slideStep.Width=2;

pard.slideStep.object=struct('Style','edit','String',5);
pard.slideStep.position=[2,3];
pard.slideStep.TooltipString = 'If you set it as 5, it means during the density comparison, every grid value will be campared to its following 4 (5 minus 1, which means the reference grid itself) right neighbors';

pard.t_adjM.object=struct('Style','text','String','Adjustment of M');
pard.t_adjM.position=[3,1];
pard.t_adjM.Width=2;

pard.adjM.object=struct('Style','edit','String',1.0003);
pard.adjM.position=[3,3];
pard.adjM.TooltipString = 'If you set it as 1.0003, it means during the optimization, if the measurment of current step (Mcur) is 0.0003-time worse than the measurment of the previous step (Mpre), Mcur will still be considered as a good result. The measurment, which defines the boundary is good or not, can be definde by users.';

% pard.dxt.Width=3;
pard.inputParameters={'numberOfLayers','sr_layerson','se_cellfov','se_sitefov','se_siteroi'};
pard.plugininfo.type='ROI_Evaluate';
end

function M = calMeasurement(x,y,qx,qy,Size)
    ref = interp1(x,y,qx);
    rightIdx = qy > ref;
    leftIdx = ~rightIdx;
    leftA = sum(y);
    rightA = Size(1)*Size(2)-sum(y);
    leftD = sum(leftIdx)/leftA;
    rightD = sum(rightIdx)/rightA;
    M = 0-((leftD-1)^2 + (rightD-0)^2)^(1/2);
end

function [newBoundary, idx] = mergeStallsByMove(boundary)
    % Merge steps by move
    move = diff(boundary(:,2));
    stalled = move==0;
    stalled = [boundary(find(stalled), 1) boundary(find(stalled)+1, 1)];
    indConStalled = find(diff(stalled,1,2)>1);
    conStalled = stalled(indConStalled,:);
    conStallTime = diff(conStalled,1,2);

    smallMove = move == 1;
    smallMove = [boundary(find(smallMove), 1) boundary(find(smallMove)+1, 1)];

    afterStall = ismember(smallMove(:,1),conStalled(:,2));
    beforeStall = ismember(smallMove(:,2), conStalled(:,1));
    abStallSM = smallMove(afterStall+beforeStall==2,:);
    lPosition2change = [];
    for i = 1:size(abStallSM,1)
    upstreamStall = conStalled(:,2) == abStallSM(i,1);
    downstreamStall = conStalled(:,1) == abStallSM(i,2);
        if conStallTime(upstreamStall) >= conStallTime(downstreamStall);
            lPosition2change = ismember(boundary(:,1), conStalled(downstreamStall,:)');
            lPositionValue = ismember(boundary(:,1), conStalled(upstreamStall,:)');
            boundary(lPosition2change,2) = boundary(lPositionValue,2);
        else
            lPosition2change = ismember(boundary(:,1), conStalled(upstreamStall,:)');
            lPositionValue = ismember(boundary(:,1), conStalled(downstreamStall,:)');
            boundary(lPosition2change,2) = boundary(lPositionValue,2);
        end
    end
    newBoundary = boundary;
    idx=lPosition2change;
end
function [newBoundary, idx] = mergeByMove(boundary)
    move = diff(boundary(:,2));
    point2drag = move == 1;
    indPoint2drag = find(point2drag)+1;
    boundary(indPoint2drag,2) = boundary(indPoint2drag-1,2);
    newBoundary = boundary;
    idx = indPoint2drag;
end
function newBoundary = mergeByTime(boundary)
    point2keepCp = [];
    while 1
        stallTime = diff(boundary(:,1));
        point2keep = stallTime > 5;
        if isequal(point2keep, point2keepCp)
            break
        end
        point2keepCp = point2keep;
        indPoint2keep=unique([find(point2keep); find(point2keep)+1]);
        point2keep(indPoint2keep)=1;
        boundary = boundary(point2keep,:);
        newBoundary = boundary;
    end
end

function salltime = calStallTime(boundary)
    move = diff(boundary(:,2));
    stalled = move==0;
    stalled = [boundary(find(stalled), 1) boundary(find(stalled)+1, 1)];
    salltime = diff(stalled,1,2);    
end

function growthTime = calGrowthTime(boundary)
    move = diff(boundary(:,2));
    growth = move>0;
    growth = [boundary(find(growth), 1) boundary(find(growth)+1, 1)];
    growthTime = diff(growth,1,2);    
end

function stepWidth = calStepWidth(boundary)
    move = diff(boundary(:,2));
    stepWidth = move(move>0);
end

function avgRate = calAvgRate(boundary)
    startNEnd = boundary([1 end], :);
    diffTimeNPosition = diff(startNEnd, 1, 1);
    avgRate = diffTimeNPosition(2)/diffTimeNPosition(1);
end