function  savesml(locData,file,p)

if p.saveroi
    [~,indgroi]=locData.getloc('xnm','position','roi');  
    indgl=false(size(indgroi));
    for k=1:p.numberOfLayers
        if p.sr_layerson(k)
           [~,indgh]=locData.getloc('xnm','layer',k); 
           if length(indgh)~=length(indgl)
               disp('save visible works only for ungrouped data')
           end
           indgl=indgl|indgh;
       
        end
    end
    indg=indgl&indgroi;
    
else
    indg=[];
end
saveloc=locData.savelocs([],indg); % BETA , maybe problematic with more than 1 file: this will save only displayed loicalizations
% if ~isempty(locData.SE)
%     saveloc.siteexplorer=locData.SE.save;
% end


if locData.getGlobalSetting('saveas73')
    version='-v7.3';
else
    version='-v7';
end
rg=p.mainGui; 
parameters=rg.saveParameters;
fileformat.name='sml';
out=struct('saveloc',saveloc,'fileformat',fileformat,'parameters',parameters);
if isfield(locData.files.file(1),'transformation')
    for k=length(obj.locData.files.file):-1:1
        if ~isempty(obj.locData.files.file(end).transformation)
            out.transformation=locData.files.file(k).transformation;
            break
        end
    end
    
end
v=saverightversion(file,out,version);
disp(['saved as version ' v])
% save(file,'saveloc','fileformat','parameters','-v7.3');
end