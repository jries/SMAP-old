classdef CME2CSide_yule<interfaces.SEEvaluationProcessor
    properties
        boundary
    end
    methods
        function obj=CME2CSide_yule(varargin)        
                obj@interfaces.SEEvaluationProcessor(varargin{:});
        end
        function out=run(obj, inp)
            locs=obj.getLocs({'locprecnm','PSFxnm','xnm','ynm','frame','channel'});  
            oriCor = [locs.xnm locs.ynm];
            
            % Do PCA together
            [U,SigmaV,lambda] = pca(oriCor);
            
            
            krot=2; % Get axies
            [Urot,T] = rotatefactors(U(:,1:krot)); % get the rotation factor
            
            % Rotate the data
            Vrot=oriCor*Urot;  
            
            minVrot = min(Vrot); % Get minmun of the dataset

            sVrot = Vrot-minVrot; % shift all points to make the points attach to the x and y axies
            
             % Get the kernel matrix
            
            [inLoc1, kf1] = kernelFeatures(sVrot, locs.channel==1);
            [inLoc2, kf2] = kernelFeatures(sVrot, locs.channel==2);
            

            
            
            figure(556)
            clf
            imagesc(z1)
            hold on
            
            
            figure(559)
            clf
            [C1,H1] = contour(z1, 200);
            hold on
            L1 = 70;
            contour(z1,[H1.LevelList(L1), H1.LevelList(L1)],'LineWidth',2 );
            
          
            
            scatter(j,i)
            scatter(sVrot(:,1), sVrot(:,2), 'MarkerFaceColor',[0 0 0])
            
            
            
            
             z2 = getKernelMatrix(sVrot(locs.channel==2,:), 'Bandwidth', bw2*0.7);
            
            Y2 = quantile(z2(:),0.95);
            z2(z2<=Y2)=0;
            
            
            
            

            
            figure(555)
            clf
            imagesc(z2)
            hold on
            BW = imregionalmax(z2);
            imagesc(BW)
            contour(z2, 200)
            contour(z1, 200)
            
            scatter(sVrot(:,1), sVrot(:,2), 'MarkerFaceColor',[0 0 0])
            
            
            
            
            % offset
            Y = mean(z(:))*1.5;
            z(z<=Y)=0;

            figure(556)
            clf
            BW = imregionalmax(z);
            hold on
            imagesc(z)
            scatter(sVrot(:,1), sVrot(:,2), 'MarkerFaceColor',[0 0 0])

        end

        function pard=guidef(obj)
            pard=guidef;
        end
    end

end


function [z, bw] = getKernelMatrix(matrixSize, sVrot, varargin) 
    [meshx, meshy] = meshgrid(0:matrixSize(1), 0:matrixSize(2)); % make full grid
    q = [meshx(:), meshy(:)]; % convert the grid to positions
    if length(varargin)==0
        [sVrotx,xy,bw]=ksdensity(sVrot, q); % get kernel density estimation
    else
        [sVrotx,xy,bw]=ksdensity(sVrot, q, varargin{1,1}, varargin{1,2}); % get kernel density estimation    
    end
    
    % converted into an image
    idx = sub2ind(matrixSize(end:-1:1)+1, xy(:,2)+1, xy(:,1)+1); 
    z = zeros(matrixSize(end:-1:1)+1);
    z(idx) = sVrotx;
end
            
function pard=guidef
pard.t_gridRange.object=struct('Style','text','String','Bin size of the grids');
pard.t_gridRange.position=[1,1];
pard.t_gridRange.Width=2;
end

function filtered = corFiltedByMask(mask, cor)
    sub = find(mask);
    ccor = ceil(cor);
    q = sub2ind(size(mask), ccor(:,2)'+1, ccor(:,1)'+1);
    filtered = ismember(q,sub);
    filtered = cor(filtered,:);
end

function maxSVrot = kernelMatrixSize(sVrot)
    maxSVrot = max(sVrot);
    maxSVrot = ceil(maxSVrot); % get maximun of x and y position
end

function [filteredCor, kf] = kernelFeatures(parentalLocs, sub)
    kf = [];
    [~, bw1] = getKernelMatrix(kernelMatrixSize(parentalLocs), parentalLocs(sub,:));

    z1Ori = getKernelMatrix(kernelMatrixSize(parentalLocs), parentalLocs(sub,:), 'Bandwidth', bw1*0.7);
     % offset
    Y1 = quantile(z1Ori(:),0.8);
    z1 = z1Ori;
    z1(z1Ori<=Y1)=0;
    filteredCor = corFiltedByMask(z1, parentalLocs(sub,:));
    [~, bw1] = getKernelMatrix(kernelMatrixSize(parentalLocs), filteredCor);
    z1 = getKernelMatrix(kernelMatrixSize(parentalLocs), filteredCor, 'Bandwidth', bw1*0.7);
    
    BW = imregionalmax(z1);
    [i,j]=find(BW>=1);

    peakVal = z1(sub2ind(size(BW), i, j));
    lessThanf = find(max(peakVal)./peakVal<=1.2);
    peakLoc = [j i];
    peakLoc = peakLoc(lessThanf,:);
    numPeak = length(lessThanf);


    Q1 = quantile(filteredCor,0.25);
    Q2 = quantile(filteredCor,0.5);
    Q3 = quantile(filteredCor,0.75);
    DQ = Q3-Q1;
    shiftFactor = (Q2-Q1)./(Q3-Q2);
    
    pd = pdist(peakLoc, 'euclidean');
    
    kf.Q1 = Q1;
    kf.Q2 = Q2;
    kf.Q3 = Q3;
    kf.DQ = DQ;
    kf.shiftFactor = shiftFactor;
    kf.numPeak = numPeak;
    kf.peakLoc = peakLoc;
    kf.peakDist = pd;
end