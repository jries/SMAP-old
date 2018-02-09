function tom_makefun
%TOM_MAKEFUN compiles C files
%
%   tom_makefun
%
%   Compiles the C files which are underlying some TOM toolbox functions.
%   Run the function on your computer system, if the C files are not
%   precompiled for your computer system. The C files are precompiled
%   for Linux, Windows 32bit and Windows 64bit.
%   If you want to run TOM_MAKEFUN it is necessary to move the MATLAB
%   current directory to your TOM folder.
%
%EXAMPLE
%   cd /fs/tom;
%   tom_makefun;
%
%NOTE
% works only in the tom folder 
%
%
%REFERENCES
%
%SEE ALSO
%   MEX
%
%   created by ME 01/02/08
%
%   Nickell et al., 'TOM software toolbox: acquisition and analysis for electron tomography',
%   Journal of Structural Biology, 149 (2005), 227-234.
%
%   Copyright (c) 2004-2007
%   TOM toolbox for Electron Tomography
%   Max-Planck-Institute of Biochemistry
%   Dept. Molecular Structural Biology
%   82152 Martinsried, Germany
%   http://www.biochem.mpg.de/tom

cd IOfun/
mex tom_emreadinc.c
disp('Function tom_emreadinc.c compiled!');
mex tom_emreadinc_resample2.c
disp('Function tom_emreadinc_resample2.c compiled!');
cd ..
cd Reconstruction/
mex tom_backproj3dc.c
disp('Function tom_backproj3dc.c compiled!');
mex tom_vol2projc.c
disp('Function tom_vol2projc.c compiled!');
cd ..
cd  Sptrans/
mex tom_rotatec.c
disp('Function tom_rotatec.c compiled!');
cd ..

