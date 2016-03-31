dirlist=genpath('shared');
addpath(dirlist)
if ~isempty(who('g'))
    delete(g)
end
 close all
 
g=gui.GuiMainSMAP;g.makeGui;    