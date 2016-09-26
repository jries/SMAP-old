function io=imageloaderAll(file)
% imageloaderAll selects the right image loader based on the presence of
% metadata.
   [path,~,ext]=fileparts(file);
   switch ext
       case '.tif'
           if exist([path filesep 'metadata.txt'],'file')
               imloader=@imageloaderMMsingle;
           elseif ~isempty(dir([path filesep '*metadata.txt']))
                imloader=@imageloaderMMstack;
           elseif ~isempty(strfind(file,'img_'))
               imloader=@imageloaderMMsingle;
           else
               imloader=@imageloaderOME;
           end
       otherwise
           imloader=@imageloaderOME;
   end    
   io=imloader(file);
end



%distinguish:
    %MM Tif single (with Metadata)
    %MM Tif stack (with Metadata)
    %Tif single wihtout metadata
    %Tif stack without metadata
    %any ome compatible file