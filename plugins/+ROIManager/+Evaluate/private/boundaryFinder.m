% Import the kimograph
global se
sites=se.sites;
test = sites(1).evaluation.BALM_fibril_growth.kimograph;
test(test>=1)=1;

% Binarize the kimograph
[testRefx, testRefy] = find(test);
[testNZx, testNZy] = find(test);

% use euclidean distances to filter out noises
ED = pdist2([testNZx testNZy],[testNZx testNZy],'euclidean');
sortedED = sort(ED,2);
sortedED = sortedED(:,1:20);
rmean = mean(sortedED,2);
test(sub2ind(size(test), testNZx([rmean>=5]), testNZy([rmean>=5]))) = 0;

% Get cordinates of the kemograph
Size = size(test);
[CorX CorY] = meshgrid(1:Size(1), 1:Size(2));
Order = sub2ind(Size, CorX(:), CorY(:));

% Make grids
gridRange = 5;
adjCorX = ceil(CorX(:)./gridRange);
adjCorY = ceil(CorY(:)./gridRange);

% Get density
D = accumarray([adjCorX adjCorY], test(Order), [], @sum);

% Filter out grids that might be background
cutoff = 3;
D(D<=cutoff)=0;

% Do sliding subtraction
SizeD = size(D);
slideStep = 5;

diff = cell(slideStep-1,1);
diffAll = zeros(SizeD(1), SizeD(2)-slideStep);
for i = 1:(slideStep-1)
    diffAll = diffAll + (D(:,1+i:end-slideStep+i) - D(:,1:(end-slideStep)) < 0);
end

% Get grids with 4 zero grids on its right-hand side
diffAllLog = diffAll == 4; 
[bX,bY] = find(diffAllLog, 1, 'last');

% Get boundary (itself is 1, and its cumsum is also 1)
boundary = diffAllLog & cumsum(diffAllLog,2,'reverse') == 1;

%figure(50)
%mesh(boundary)
%rotate3d on

% Find the corrdinates of boundary
[Cx, Cy]= find(boundary,sum(boundary(:)),'last');

% Make sure every frame (x-axies) have a point
CxFinal = (1:75)';
CyFinal = zeros(75,1);
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




figure(50)
scatter(testNZx, testNZy)
hold on
plot(CxFParent, CyFParent)


Mref = calMeasurement(CxFParent,CyFParent,testNZx,testNZy,Size);
MrefO = Mref;
CyFParentA = CyFParent;
for i=size(CxFParent):-1:1
   thisFrameY = testNZy(testNZx == CxFParent(i));
   thisFrameY = thisFrameY(thisFrameY>CyFParent(i));
   ii = 1;
   while ii <= length(thisFrameY)
       CyFParentO = CyFParentA;
       
       % New boundary
       CyFParentA(i)=thisFrameY(ii);
       
       CyFParentA = cummax(CyFParentA);

       % fit to new boundary
       Mnew = calMeasurement(CxFParent,CyFParentA,testNZx,testNZy,Size);
       if Mnew >= MrefO*1.0003 && Mnew >= Mref
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
figure(51)
scatter(testNZx, testNZy)
hold on
plot(CxFParent, CyFParentA)
hold off
end

function M = calMeasurement1(x,y,qx,qy,Size)
    ref = interp1(x,y,qx);
    rightIdx = qy > ref;
    leftIdx = ~rightIdx;
    B = sum(y);
    M = sum(leftIdx)/B;
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

function M = calMeasurement0(x,y,qx,qy,Size)
    ref = interp1(x,y,qx);
    rightIdx = qy > ref;
    leftIdx = ~rightIdx;
    leftA = polyarea([qx(end);x], [0; y]);
    rightA = polyarea([x;0], [y; qy(end)]);

    M = (sum(leftIdx)/leftA)/(sum(rightIdx)/rightA);
end