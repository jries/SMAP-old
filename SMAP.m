% t=timerfindall;
% delete(t);

dirlist=genpath('shared');
addpath(dirlist)
% rootl=fileparts(pwd);
% addpath([rootl filesep 'mYlibraries' filesep 'myfunctions']);
% addpath([rootl filesep 'mYlibraries' filesep 'myMATLABfunctions']);
% addpath([rootl filesep 'mYlibraries' filesep 'external' filesep 'externaltools']);
% addpath([rootl filesep 'mYlibraries' filesep 'external' filesep 'export_fig']);
if ~isempty(who('g'))
    delete(g)
end
 close all
 
g=gui.GuiMainSMAP;g.makeGui;    