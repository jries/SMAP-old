dirlist=genpath('shared');
addpath(dirlist)
if ~isempty(who('g'))
    delete(g)
end
%  close all
try
g=gui.GuiMainSMAP;g.makeGui;       
catch err
    disp('Error making the GUI. Try deleting plugins/plugin.m and the settings/temp directory.')
 
%     g.status(err.message)
    err.rethrow
end

[status,message]=system('git status');
if status==0
    ind=find(message==10);
    disp(message(1:ind(2)));
end

 