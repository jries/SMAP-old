function i = randindex(objs)
% I = RANDINDEX(OBJS) returns a random index into a vector or cell
% array OBJS of objects.  If OBJS is a matrix, RANDINDEX returns a random
% row index from the matrix.

if size(objs, 1) > 1
    i = 1+fix(rand*size(objs,1));
else
    i = 1+fix(rand*length(objs));
end

