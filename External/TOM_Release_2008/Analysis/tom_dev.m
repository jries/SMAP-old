function [a,b,c,d,e]=tom_dev(A,suppress)
%TOM_DEV Calculates the mean, max, min, standard-deviation, variance
%   of an input image
%
%   [a,b,c,d,e]=tom_dev(A,suppress)
%
%PARAMETERS
%  INPUT
%   a           ...
%   suppress    ...
%   IN      : input image
%   flag    : This is optional.
%               -'noinfo' Not display the information
%
%   mean    : Mean value
%   max     : maximum value
%   min     : minimum value
%   Std     : standard deviation
%   variance: Variance of image
%
%  OUTPUT
%   a           ...
%   b           ...
%   c           ...
%   d           ...
%   e           ...
%
% EXAMPLE
%   [mean, max, min, std, variance] = tom_dev(IN,flag)
%   im = tom_emread('proteasome.em');
%   tom_dev(im.Value);
%
%REFERENCES
%
%SEE ALSO
%   ...
%
%   created by AF 03/18/02
%   updated by WDN 07/11/03
%   updated by 05/23/06
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

%error(nargchk(1,2,nargin))
[s1,s2,s3,s4]=size(A);
a=sum(sum(sum(sum(A))))./(s1.*s2.*s3.*s4);
b=max(max(max(max(A))));
c=min(min(min(min(A))));
d=std2(A);
e=d.^2;
if nargin <2
    f=sprintf('Mean= %g,  Max= %g,  Min= %g,  Std= %g,  Variance= %g', a,b,c,d,e);disp(f);
elseif nargin==2
    switch suppress
        case 'noinfo'
            %no action
        otherwise
            f=sprintf('Mean= %g,  Max= %g,  Min= %g,  Std= %g,  Variance= %g', a,b,c,d,e);disp(f);
    end
end
