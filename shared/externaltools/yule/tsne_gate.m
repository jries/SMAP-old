function [ gated_indices ] = tsne_gate(plothandle, mapped_data, ~, isPoly)

    if (isPoly) 
%         [pl,xs,ys] = selectdata('sel','lasso');
%         gated_indices = find(mapped_data==pl);
        h = impoly(gca);
        disp('waiting: Double click on node when finished');
        vert = wait(h);
        disp('before ray casting');
        gated_indices = find(inpoly(mapped_data, vert));

%     ray casting - PIP in C++ for an individual point - project point to y-axis and
%     count an odd number of vertices is to the right of point. 
%         gated_indices = []
% 
%         % iterate over all points
%         for pInd = 1:size(mapped_data, 1)
%             p = mapped_data(pInd, :);
%             c = 0;
%             j=size(vert, 1); % j starts at the index of the last vertice
% 
%             % iterate over vertices
%             for i = 1:size(vert, 1)
%                 if ((vert(i,2)>p(2)) ~= (vert(j, 2)>p(2))) & ...
%                     (p(1) < (((vert(j,1)-vert(i,1)) * (p(2)-vert(i,2)) / (vert(j, 2)-vert(i, 2))) + vert(i, 1)))
%                     c = ~c;
%                 end
%                 j = i; % j is always one behind 1 except for the first run.
%             end
% 
%             % if point is to the left of an odd number of vertices
%             if c==1
%                 gated_indices(end+1) = pInd;
%             end
% 
% 
%         end
    else 
    
        rect = getrect(plothandle);
        left = rect(1);
        bottom = rect(2);
        width = rect(3);
        height = rect(4);

        gated_indices = find((mapped_data(:,1)>left) & (mapped_data(:,1)<left+width) & (mapped_data(:,2)>bottom) & (mapped_data(:,2)<bottom+height));
    end           
end

