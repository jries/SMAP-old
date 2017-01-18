function transformation=loadtransformation(obj,file,dataset)
if nargin<3
    dataset=length(obj.locData.files.file);
end
if exist(file,'file')
    l=load(file,'transformation');
    transformation=l.transformation;
elseif isfield(obj.locData.files.file(1),'transformation')
     if ~isempty(obj.locData.files.file(dataset).transformation)
            transformation=locData.files.file(dataset).transformation;
     else
        for k=length(obj.locData.files.file):-1:1
            if ~isempty(obj.locData.files.file(k).transformation)
                transformation=locData.files.file(k).transformation;
                break
            end
        end
     end  
else
    transformation=[];
end