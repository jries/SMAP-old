%%%%% The FOCAL analysis perform clustering of super resolution localization
%%%%% data. The input localization table should be in the same format of 
%%%%% RapidStorm localization  tables.
%%%%% The program has 3 modules: 
%%%%%        1. ConfineLocTable: removes localizations that fall out of the
%%%%%        image boundries. This can happen because of drift correction.
%%%%%        2. DensityFilter: only keeps localizations that are in
%%%%%        those pixels of density map that have values >=minL (core pixels) 
%%%%%        or have a core pixel in their 4cc neighbourhood (border pixels).
%%%%%        3. LocToROI: Identifies connected core and border pixels as a
%%%%%        cluster candidates. Only Those with span>=10 pixels will be
%%%%%        accepted as a cluster.
%%%%%
%%%%% The program saves 3 localization files. One after each module. 
%%%%% The 3rd file contain the list of localizations within clusters found by
%%%%% FOCAL analysis. This table can be reconstructed by rapidstorm to
%%%%% visulaize the final result. 
%%%%% 
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

pixel=100;  % size of each camera pixel in the sample plane nm
SRpixel=10; % Should set equal to localization uncertainty but ususally 10nm is an appropriate choice
minL=14;    
minC=10;   % minC= minimum size of acceptable clusters: in current version of FOCAL minC=10. 
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
% %%%%%%%  2. Cluster the core and border localizations and save the
% %%%%%%%     Localizations within accepted clusters in a file
% %%%%%%%     at [FOCALPathName FOCALFileName]
% %%%%%%%  
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      
        [DFPathName,DFFileName] = DensityFilter(BFileName,BPathName,minL);  
        [FOCALPathName,FOCALFileName]= LocToROI(DFFileName,DFPathName,minC); 
   





        