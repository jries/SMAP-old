function out = tom_isspiderfile(filename)
%TOM_ISSPIDERFILE checks for spiderfiles
%
%   out = tom_isspiderfile(filename)
%
%PARAMETERS
%
%  INPUT
%   filename            filename
%  
%  OUTPUT
%   out                 1 for spiderformat
%                       0 no spiderformat 
%
%EXAMPLE
%  
% im=tom_emread('pyrodictium_14.em');
% tom_spiderwrite('test.spi',im.Value);
% out=tom_isspiderfile('test.mrc');
%
%
%
%REFERENCES
%
%SEE ALSO
%   TOM_SPIDERWRITE, TOM_SPIDERHEADER, TOM_SPIDERREAD
%
%   created by AK 04/25/06
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


fid = fopen(filename,'rb');
if fid==-1
    out=-1;
    return;
end

test = fread(fid,5,'float');

if (isempty(test)==1)
    out=-1;
    fclose(fid);
    return;
end;
    
if test(1) > 0 && test(1) < 10000 && test(1) > 0 && test(2) < 10000 && (test(5) == 1 || test(5) == 3 || test(5) == -11 || test(5) == -12 || test(5) == -21 || test(5) == -22)
   out = 1; 
else
   out = 0;
end

fclose(fid);


