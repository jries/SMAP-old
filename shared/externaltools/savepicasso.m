function savepicasso(PathName,FileName, locs)

%savepicasso: save molecule list as hdf5 file and additional yaml file to
%be readable with picasso
%HDF5 import modified from the 'h5ex_t_cmpd.m' example from hdfgroup.com

%Input: PathName: old file path e.g. C:/users/
%Input: FileName: old FileName with ending e.g. locs.mat
%Input: struct locs:
%locs.frame:    frames, picasso frames start with 0
%locs.x:        x-coordinates in pixels, make sure that x>0
%locs.y:        y-coordinates in pixels, make sure that y>0
%locs.lpx:      localization precision in x in pixels
%locs.lpy:      localization precision in y in pixels

%2017 m.strauss@biochem.mpg.de

%Make sure all fields are double 

locs.frame = double(locs.frame);
locs.x = double(locs.x);
locs.y = double(locs.y);
locs.lpx = double(locs.lpx);
locs.lpy = double(locs.lpy);

%get parameters for yaml
width = ceil(max(locs.x));
height = ceil(max(locs.y));
frames = max(locs.frame);

%Create and hdf5 file with /locs
saveName       = fullfile(PathName,strcat(FileName,'.mat.hdf5'));
DATASET        = '/locs';
dims          = length(locs.frame);


file = H5F.create (saveName, 'H5F_ACC_TRUNC',...
    'H5P_DEFAULT', 'H5P_DEFAULT');

% Assign the required data types (all as double)

doubleType=H5T.copy('H5T_NATIVE_DOUBLE');
sz(1)     =H5T.get_size(doubleType);
doubleType=H5T.copy('H5T_NATIVE_DOUBLE');
sz(2)     =H5T.get_size(doubleType);
doubleType=H5T.copy('H5T_NATIVE_DOUBLE');
sz(3)     =H5T.get_size(doubleType);
doubleType=H5T.copy('H5T_NATIVE_DOUBLE');
sz(4)     =H5T.get_size(doubleType);
doubleType=H5T.copy('H5T_NATIVE_DOUBLE');
sz(5)     =H5T.get_size(doubleType);

% Computer the offsets to each field. The first offset is always zero.
offset(1)=0;
offset(2:5)=cumsum(sz(1:4));

% Create the compound datatype for memory.
memtype = H5T.create ('H5T_COMPOUND', sum(sz));
H5T.insert (memtype,...
    'frame',offset(1), doubleType);
H5T.insert (memtype,...
    'x',offset(2), doubleType);
H5T.insert (memtype,...
    'y',offset(3), doubleType);
H5T.insert (memtype,...
    'lpx',offset(4), doubleType);
H5T.insert (memtype,...
    'lpy',offset(5), doubleType);

%Create the compound datatype with the required fields - 'frame', 'x', 'y',
%'lpx' %'lpy'
filetype = H5T.create ('H5T_COMPOUND', sum(sz));
H5T.insert (filetype, 'frame', offset(1),doubleType);
H5T.insert (filetype, 'x', offset(2), doubleType);
H5T.insert (filetype, 'y',offset(3), doubleType);
H5T.insert (filetype, 'lpx',offset(4), doubleType);
H5T.insert (filetype, 'lpy',offset(5), doubleType);

% Create dataspace.  Setting maximum size to [] sets the maximum
% size to be the current size.
space = H5S.create_simple (1,fliplr( dims), []);

% Create the dataset and write the compound data to it.
dset = H5D.create (file, DATASET, filetype, space, 'H5P_DEFAULT');
H5D.write (dset, memtype, 'H5S_ALL', 'H5S_ALL', 'H5P_DEFAULT', locs);

% Close and release resources.
H5D.close (dset);
H5S.close (space);
H5T.close (filetype);
H5F.close (file);

%Also create a yaml file for picasso to import
yamlName       = fullfile(PathName,strcat(FileName,'.mat.yaml'));
fileID = fopen(yamlName,'wt');
fprintf(fileID,'Byte Order: <\n'); %the byte-order  and data type should not matter at this point as the raw data is already localizaed
fprintf(fileID,'Data Type: uint16 \n');
fprintf(fileID,['Frames: ',num2str(frames),' \n']);
fprintf(fileID,['Width: ',num2str(width),' \n']);
fprintf(fileID,['Height: ',num2str(height),' \n']);
fprintf(fileID,'Generated by: Picasso Convert\n');
fclose(fileID);
end