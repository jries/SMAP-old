function y = scale(ival, fval, frac)
% Y = SCALE(IVAL, FVAL, FRAC) temporal scaling function
y = ival + frac * (fval - ival);
