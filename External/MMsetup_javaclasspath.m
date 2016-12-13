function [] = MMsetup_javaclasspath(path2MM)
fileList = getAllFiles(path2MM);
fileListJarBool = regexp(fileList,'.jar$','end');
fileListJarBool = cellfun(@isempty,fileListJarBool);
fileListJar = fileList(~fileListJarBool);
fid = fopen(fullfile(prefdir,'MMjavaclasspath.txt'),'w');
fprintf(fid,'<before>\r\n');
cellfun(@(x) fprintf(fid,'%s\r\n',x), fileListJar);
fclose(fid);
%% nested directory listing ala gnovice from stackoverflow
% inputs and outputs are self-explanatory
function fileList = getAllFiles(dirName)
dirData = dir(dirName);      % Get the data for the current directory
dirIndex = [dirData.isdir];  % Find the index for directories
fileList = {dirData(~dirIndex).name}';  % Get a list of the files
if ~isempty(fileList)
    fileList = cellfun(@(x) fullfile(dirName,x),fileList,'UniformOutput',false);
end
subDirs = {dirData(dirIndex).name};  % Get a list of the subdirectories
validIndex = ~ismember(subDirs,{'.','..'});  % Find index of subdirectories
%   that are not '.' or '..'
for iDir = find(validIndex)                  % Loop over valid subdirectories
    nextDir = fullfile(dirName,subDirs{iDir});    % Get the subdirectory path
    fileList = vertcat(fileList, getAllFiles(nextDir));  % Recursively call getAllFiles
end
