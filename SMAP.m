dirlist=genpath('shared');
addpath(dirlist)
if ~isempty(who('g'))
    delete(g)
end
%  close all
try
g=gui.GuiMainSMAP;g.makeGui;    
catch err
    disp('try deleting plugins.m in the plugins directory and the temp directory in the settings directory')
    err.rethrow
end

[status,message]=system('git status');
if status==0
    ind=find(message==10);
    disp(message(1:ind(2)));
end

