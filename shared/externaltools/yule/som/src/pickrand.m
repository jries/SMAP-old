function obj = pickrand(objs)
% OBJ = PICKRAND(OBJS) returns a random object OBJ from a vector or cell
% array OBJS of objects.  If OBJS is a matrix, PICKRAND returns a random
% row from the matrix.

i = randindex(objs);

if size(objs, 1) > 1
    obj = objs(i,:);
else
    obj = objs(i);
end

