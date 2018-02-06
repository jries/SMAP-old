% Import the kimograph
test = sites(2).evaluation.BALM_fibril_growth.kimograph;

% Get cordinates of the kemograph
Size = size(test);
[CorX CorY] = meshgrid(1:Size(1), 1:Size(2));
Order = sub2ind(Size, CorX(:), CorY(:));

% Make grids
gridRange = 4;
adjCorX = ceil(CorX(:)./gridRange);
adjCorY = ceil(CorY(:)./gridRange);

% Get density
D = accumarray([adjCorX adjCorY], test(Order), [], @sum);

% Filter out grids that might be background
cutoff = 3;
D(D<=cutoff)=0;

% Do sliding subtraction
SizeD = size(D);
slideStep = 7;

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

figure(50)
mesh(boundary)
rotate3d on

% Find the corrdinates of boundary
[Cx, Cy]= find(boundary,sum(boundary(:)),'last');

% Make sure every frame (x-axies) have a point
CxFinal = (1:75)';
CyFinal = zeros(75,1);
CyFinal(Cx) = Cy;

% Find the cummax (next points should always >= the previous points)
CyFinal = cummax(CyFinal);
CxFParent = CxFinal*gridRange;
CyFParent = CyFinal*gridRange;

[testNZx, testNZy] = find(test);

figure(50)
scatter(testNZx, testNZy)
hold on
plot(CxFParent, CyFParent)