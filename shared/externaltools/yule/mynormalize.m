function data = mynormalize(data, percentile, varargin)
% normalize data according to specified percentile
    
    for i=1:length(varargin)-1
        if(strcmp(varargin{i},'num_locs'))
            num_locs = varargin{i+1};
        end
    end
    
    fprintf('Normalizing according to the %gth percentile... \n', percentile);
    data = data-repmat(prctile(data, 100-percentile, 1), size(data,1),1);
    data = data./repmat(prctile((data), percentile, 1),size(data,1),1);

    data(data > 1) = 1;
    data(data < 0) = 0;
    data(isinf(data)) = 0;
end