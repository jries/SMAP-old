function plot_points(a)
% PLOT_POINTS(A) displays points of Nx2 or Nx3 matrix A in the unit square
% or cube.

if size(a, 2) == 3
    plot3(a(:,1), a(:,2), a(:,3), '.', 'MarkerSize', 1)
    axis([0 1 0 1 0 1])

else
    plot(a(:,1), a(:,2), '.', 'MarkerSize', 1)
    axis([0 1 0 1])
end

axis square
