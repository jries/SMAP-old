function out=tom_isxmippsell(filename)
%TOM_ISXMIPPSELL checks for xmippsell
%
%   out=tom_isxmippsell(filename)
%
%PARAMETERS
%
%  INPUT
%   filename            name of the cell-file
%  
%  OUTPUT
%   out                 0/1 
%
%EXAMPLE
%   
%      
%
%REFERENCES
%
%SEE ALSO
%   tom_imagicread, tom_imagicwrite , tom_isspiderfile, tom_isemfile
%
%   created by fb 28/01/06
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


flag=0;

[path name ext]=fileparts(filename);

if (strcmp(ext,'.sel')==0)
    out=0;
    return;
end;


fid=fopen(filename);
line=fgetl(fid);
fclose(fid);

line=strtok(line,' ');

if (fid==-1)
    out=0;
    return;
end;

try 
    im=tom_spiderread(line);
    flag=1;
catch
    disp([line ' not found. Change directory.']);
end;
 

try
    im=tom_spiderread(['.' line]);
    flag=2;    
catch

end;

if (flag==0)
    out=0;
    return;
end;


out=1;


