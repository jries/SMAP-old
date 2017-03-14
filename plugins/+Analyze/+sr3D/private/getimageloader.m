function il=getimageloader(obj,filename)
filename=strrep(filename,'\','/');
[~,~,ext]=fileparts(filename);

if isempty(ext)&&~strcmp(filename(end),'/') %directory
    filename=[filename '/'];
end
try
    il=imageloaderAll(filename,[],obj.P); %still exist?
catch
    maindir=obj.getGlobalSetting('DataDirectory');
    filenamef=findfilepath(filename,maindir);
    try
        il=imageloaderAll(filenamef,[],obj.P); %still exist?
    catch
        lastsml=obj.getPar('lastSMLFile');
        filenamef=findfilepath(filename,fileparts(lastsml));
            try
                il=imageloaderAll(filenamef,[],obj.P); %still exist?
            catch
                [~,~,ext]=fileparts(filename);
                if isempty(ext)
                    filenamef=[filename '.tif'];
                end
                [f,path]=uigetfile(filenamef);
                if f
                    il=imageloaderAll([path f],[],obj.P); 
                else
                    il=[];
                end
            end
    end
    %look for file in main directory
end
end

function filename=findfilepath(filename,maindir)
 filename=strrep(filename,'\','/');
    d=dir(maindir);
    alldir={d([d.isdir]).name};
    ind=[1 strfind(filename,'/')];
    for k=1:length(ind)-1
        thisf=filename(ind(k)+1:ind(k+1)-1);
        if any(strcmp(alldir,thisf))&&~isempty(thisf)


            filename=[maindir '/' filename(ind(k)+1:end)];
            break
        end
    end

end