function [G1, S, G2, M0, M1] = cyclerclassification_v2(featurenames, data)

debug = false;

%% extract 5 features
ind = cellstrfnd(featurenames, 'Blue') &...
    cellstrfnd(featurenames, 'Nuclei') &...
      cellstrfnd(featurenames, 'IntegratedIntensity') &...
      ~cellstrfnd(featurenames, 'Edge');
dapi = data(:, ind);

ind = cellstrfnd(featurenames, 'Red') &...
    cellstrfnd(featurenames, 'Nuclei') &...
      cellstrfnd(featurenames, 'IntegratedIntensity') &...
      ~cellstrfnd(featurenames, 'Edge');
edu = data(:, ind);

ind = cellstrfnd(featurenames, 'Blue') &...
      cellstrfnd(featurenames, 'Texture') &...
      cellstrfnd(featurenames, 'Nuclei') &...
      cellstrfnd(featurenames, 'Variance')&...
      ~cellstrfnd(featurenames, 'SumVariance')&...
      ~cellstrfnd(featurenames, 'DifferenceVariance');
dapi_text1 = data(:, ind);

ind = cellstrfnd(featurenames, 'Blue') &...
      cellstrfnd(featurenames, 'Texture') &...
      cellstrfnd(featurenames, 'Nuclei') &...
      cellstrfnd(featurenames, 'SumVariance')&...
      ~cellstrfnd(featurenames, 'DifferenceVariance');
dapi_text2 = data(:, ind);

ind = cellstrfnd(featurenames, 'Red') &...
      cellstrfnd(featurenames, 'Texture') &...
      cellstrfnd(featurenames, 'Nuclei') &...
      cellstrfnd(featurenames, 'InfoMeas1');
edu_text = data(:, ind);

cyclerfeatures = [dapi edu dapi_text1 dapi_text2 edu_text];

if size(cyclerfeatures, 2) ~= 5
    % we are missing one of the needed features by name
	msgID = 'cyclerclassification_v2:ccfeatures';
    msg = 'deafult 5 cell cycle features were not found by name.';
    causeException = MException(msgID,msg);   
    throw(causeException);
end

%% run emgm on 5 features

% init the emgm clustering
ninit = mynormalize([dapi edu dapi_text1], 98);
Si  = knnsearch(ninit, [.5, .75 0]);
g1i = knnsearch(ninit, [.25 .25 0]);
g2i = knnsearch(ninit, [.75 .25 0]);
mi  = knnsearch(ninit, [.5 .25 1]);

[cellcycleemgm, model, llh] = emgm((cyclerfeatures)', cyclerfeatures([Si g1i g2i mi], :)');
cellcycleemgm = cellcycleemgm';

if debug
    % sanity check
    cellcycle = data(:, find(cellstrfnd(featurenames, 'CellCycle')));
    well      = data(1, find(cellstrfnd(featurenames, 'wellid')));
    
    dprc = prctile(dapi, [2,98]);
    eprc = prctile(edu, [2,98]);
    c=distinguishable_colors(8);

    f=figure(20);
    subplot('position', [0.05 0.05 .4 .9]);
    scatter(dapi, edu, 100, cellcycleemgm, '.');
    colormap(c(1:4,:));
    axis([dprc eprc]);
    ri = clusteringError(cellcycle,cellcycleemgm)
    title(sprintf('RI:%1.2f', ri));
    subplot('position', [.55 0.05 .4 .9]);
    scatter(dapi, edu, 100,cellcycle+4,'.');
    colormap(c(5:8,:));
    axis([dprc eprc]);
    drawnow;
    screen2eps(f, sprintf('cellclassificationbench-well-%g.eps', well));
end

%% assign cell cycle label to each cluster
mean_edu = accumarray(cellcycleemgm, edu)./accumarray(cellcycleemgm, ones(numel(cellcycleemgm),1));
mean_dapi = accumarray(cellcycleemgm, dapi)./accumarray(cellcycleemgm, ones(numel(cellcycleemgm),1));
mean_dapitext = accumarray(cellcycleemgm, dapi_text1)./accumarray(cellcycleemgm, ones(numel(cellcycleemgm),1));
[~, Si] = max(mean_edu);
[~, Mi] = max(mean_dapitext);
G = setdiff(1:4,[Mi Si]);
[~, G1i] = min(mean_dapi(G));
[~, G2i] = max(mean_dapi(setdiff(1:4,[Mi Si])));
G1i = G(G1i);
G2i = G(G2i);

%% find low density cavity between G1 and G2
dapic = dapi(inrange(dapi, [prctile(dapi, 5),prctile(dapi, 95)]) & inrange(edu, [prctile(edu, 5),prctile(edu, 95)]));
educ = edu(inrange(dapi, [prctile(dapi, 5),prctile(dapi, 95)]) & inrange(edu, [prctile(edu, 5),prctile(edu, 95)]));
dapieduc = [dapic educ];

band = 32;
[~,density,X,Y] = kde2d(dapieduc, band, min(dapieduc), max(dapieduc));
cavities = density<prctile(density(:), 10);
if debug
    scatter(X(cavities), Y(cavities), 150, '.r');
end

% lower continous region for removal
[conts,regions] = bwboundaries(cavities, 'noholes');

coms = zeros(numel(conts),2);
for i=1:numel(conts)
    if (size(conts{i},1) > 5)
        coms(i,:) = sum(conts{i})./size(conts{i},1);
    end
end

cavityi =  knnsearch(coms, [0, size(regions, 2)/2]);

% remove S and cavity and use emgm+texture features to segment M\G1\G2
idx=knnsearch([X(:), Y(:)], [dapi edu]);
purge = ismember(idx, find(regions(:) == cavityi));

if debug
    scatter(dapi(purge), edu(purge), 150,'.m');
end

%% remove cavity and outliers from G1,G2,and S
S  = (cellcycleemgm==Si)  & ~purge;
G1 = (cellcycleemgm==G1i) & ~purge;
G2 = (cellcycleemgm==G2i) & ~purge;
M = (cellcycleemgm==Mi)  & ~purge;

%% filter right and left on cavity
cavityidx=knnsearch([dapi(purge) edu(purge)], [dapi edu]);
dapi_purge = dapi(purge);
divide = dapi_purge(cavityidx);
if isempty(cavityidx)
    divide = mean(prctile(dapi, [5,95]));
end

G1 = G1 & dapi < divide;
G2 = G2 & dapi > divide;
M0 = M & dapi  < divide;
M1 = M & dapi  > divide;

%% filter G1 G2 and S for outliers on dapi\edu
G1 = find(G1);
G1 = G1(~debris(cyclerfeatures(G1, 1:2)));

S = find(S);
S = S(~debris(cyclerfeatures(S, [1 2 5])));

G2 = find(G2);
G2 = G2(~debris(cyclerfeatures(G2, 1:2)));

G1 = G1(dense(dapi(G1), edu(G1), 2.5));
S  = S(dense(dapi(S), edu(S), 10));
G2 = G2(dense(dapi(G2), edu(G2), 10));

if debug
    plotCyclePhases(dapi, edu, G1,G2,S);
end

M0 = find(M0);
M1 = find(M1);

end

function plotCyclePhases(dapi, edu, G1,G2,S)
    f=figure(5);
    scatter(dapi(G1), edu(G1), '.b');
    hold on;
    scatter(dapi(S), edu(S), '.g');
    scatter(dapi(G2), edu(G2), '.r');
end

function densei=dense(x,y, prc)
[~,density,X,Y] = kde2d([x,y], 32, min([x,y]), max([x,y]));
densityidx=knnsearch([X(:), Y(:)], [x,y]);
density = density(densityidx); % density value per cell
thresh = prctile(density, prc);
densei = (density>thresh);
end
