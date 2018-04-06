function veci=inrange(A, range)
% check if vector or matrix
[rows, cols] = size(A);

if any([rows, cols] == 1) 
    veci = A>range(1) & A<range(2);
else
    if (cols ~= size(range,2))
        error('number of columns in A should be equal to number of ranges');
    end
    
    veci = zeros(rows, cols);
    for i=1:cols
        veci(:,i) = inrange(A(:,i), range(:,i));
    end
end
end