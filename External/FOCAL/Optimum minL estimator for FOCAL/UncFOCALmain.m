%%%%% This program performs the localization uncertainty analysis of the 
%%%%% super resolution localization data to obtain the optimum minL
%%%%% (minL*).
%%%%% The input localization table should be in the same format of 
%%%%% RapidStorm localization  tables.
%%%%% Copyright (C) 2015  By Amir Mazouchi

%%%%% This program is free software; you can redistribute it and/or
%%%%% modify it under the terms of the GNU General Public License
%%%%% as published by the Free Software Foundation; either version 2
%%%%% of the License, or (at your option) any later version.

%%%%% This program is distributed in the hope that it will be useful,
%%%%% but WITHOUT ANY WARRANTY; without even the implied warranty of
%%%%% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%%%%% GNU General Public License for more details.


clc;
clear all;

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%  Set input Values
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
global pixel SRpixel Nneighbor XSRp YSRp;

pixel=100;  %size of each camera pixel in the sample plane nm
SRpixel=10; % Should set equal to localization uncertainty but ususally 10nm is an appropriate choice
VminL=[2:30]; % range of minL to search for optimum minL(minL*)
VminC=10;%  minC= minimum size of acceptable clusters: in current version of FOCAL minC=10. You can set it as a vector such as 1:15 to built minL* vs minC curve
Nneighbor=4;  % definition of  neighbors: 4 or 8-conncted nearest neighbors. In the current version of FOCAL it is always 4. If it is not 4 or 8 the program will set it to 8

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%  Get the localization table in the rapidstorm format, find the
% %%%%%%%  SRL Image size and reject localizations out of the SRL image
% %%%%%%%  boundaries. A localization can fall outside of the boundaries of
% %%%%%%%  an image in the process of drift correction.
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[FileName,PathName,FilterIndex] = uigetfile([cd '\*.txt']);
[ XSRp, YSRp ] = SRLImageSize( FileName,PathName );   
[BPathName, BFileName] = ConfineLocTable( FileName,PathName );

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%  1. Perform a density filtering to only keep localizations in the
% %%%%%%%     core and border pixels
% %%%%%%%  2. Cluster the core and border localizations
% %%%%%%%  3. Calculate localization Uncertainty curves for 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for i=1:length(VminL)
    minL=VminL(i);                                                                                         
    fprintf(['minL = ' num2str(minL) '\n']);    
    fprintf('    minC = ');
    for j=1:length(VminC)
        minC=VminC(j);
        fprintf([num2str(minC) '  ']);
        [FOCALPathName,FOCALFileName]= LocToROI(BFileName,BPathName,minC); % Here is to reject localizations in domains smaller than minC from analysis
            approvalList=PixEQminL(FOCALFileName,FOCALPathName,minL); % Here is to obtain the index list of localizations that are in pixels (of density map) with value==minL
        [LocUnc,PCI]=UncertaintyAnalyzer(FOCALFileName,FOCALPathName,approvalList);        
        MLocUnc(i,j)=LocUnc; % Experimental localization uncertainty
        MPConfInt(i,j)=PCI; % Half of confidence interval of the fit to the average localization uncertaintiy       
    end  
    fprintf(['\n']);
end

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%  Fit the localization uncertainty curves and plot them
% %%%%%%%  You may need to 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
save([PathName '/UncertaintyCurves.mat'],'VminL', 'VminC', 'MLocUnc', 'MPConfInt');
LocUncFit(MLocUnc,MPConfInt,VminL,VminC);



        