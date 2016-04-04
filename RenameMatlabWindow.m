% change title of Matlab window to facilitate use of multiple instances
% MM 4 apr 16

jDesktop = com.mathworks.mde.desk.MLDesktop.getInstance;
jDesktop.getMainFrame.setTitle('Vrp1');
clear jDesktop;