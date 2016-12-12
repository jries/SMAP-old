function io=imageloaderAll(varargin)
% imageloaderAll selects the right image loader based on the presence of
% metadata.
file=varargin{1};
   [path,~,ext]=fileparts(file);
%    info=imfinfo(file);Tiff
   switch ext
       case '.tif'
           if exist([path filesep 'metadata.txt'],'file')
%                imloader=@imageloaderMM;
               if countfiles(file)>1 && ~(any(strfind(file,'MMStack'))||any(strfind(file,'.ome.')))
                   imloader=@imageloaderMMsingle;
               else
                   imloader=@imageloaderMM;
               end
           elseif ~isempty(dir([path filesep '*metadata.txt']))
                imloader=@imageloaderMM;
           elseif any(strfind(file,'MMStack'))
                imloader=@imageloaderMM;
           elseif countfiles(file)>1000
               imloader=@imageloaderMMsingle;
           else
               imloader=@imageloaderOME;
           end
       otherwise
           imloader=@imageloaderOME;
   end    
   io=imloader(varargin{:});
end

function numf=countfiles(file)
files=myfastdir(fileparts(file),'*.tif');
sstr=regexprep(files{1},'[0-9]*','[0-9]*');
isfile=regexp(files,sstr);
numf=sum(cell2mat(isfile));

end

%distinguish:
    %MM Tif single (with Metadata)
    %MM Tif stack (with Metadata)
    %Tif single wihtout metadata
    %Tif stack without metadata
    %any ome compatible file