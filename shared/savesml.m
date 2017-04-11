function  savesml(locData,file,p,excludesavefields)

if p.saveroi
    [~,indgroi]=locData.getloc('xnm','position','roi');  
    indgl=false(size(indgroi));
    for k=1:p.numberOfLayers
        if p.sr_layerson(k)
           [locs,indgh]=locData.getloc('inlayeru','layer',k); 
%            if length(indgh)~=length(indgl)
%                disp('save visible works only for ungrouped data')
%            end
           indgl=indgl|locs.inlayeru{1};
       
        end
    end
    indg=indgl&indgroi;
    
else
    indg=true(size(locData.loc.xnm));
end

if p.savefile
    filenumber=p.dataselect.Value;
    indg=indg&locData.loc.filenumber==filenumber;
end
if all(indg)%save all
    indg=[];
end
locData.loc=rmfield(locData.loc,excludesavefields);
saveloc=locData.savelocs([],indg); % BETA , maybe problematic with more than 1 file: this will save only displayed loicalizations
% if ~isempty(locData.SE)
%     saveloc.siteexplorer=locData.SE.save;
% end
if p.savefile
    saveloc.file=saveloc.file(filenumber);
    saveloc.history=saveloc.history(filenumber);
    saveloc.loc.filenumber=ones(size(locData.loc.filenumber));
end

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
    for k=length(locData.files.file):-1:1
        if ~isempty(locData.files.file(k).transformation)
            out.transformation=locData.files.file(k).transformation;
            break
        end
    end
    
end
v=saverightversion(file,out,version);
disp(['saved as version ' v])
% save(file,'saveloc','fileformat','parameters','-v7.3');
end