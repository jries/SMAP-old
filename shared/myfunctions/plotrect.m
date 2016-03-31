function plotrect(axis,pos,att)
%pos: diagonal corners
line([pos(1) pos(3) pos(3) pos(1) pos(1)],[pos(2) pos(2) pos(4) pos(4) pos(2)],'Parent',axis,'Color',att)
% plot([pos(1) pos(3) pos(3) pos(1) pos(1)],[pos(2) pos(2) pos(4) pos(4) pos(2)],att,'Parent',axis)