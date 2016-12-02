function [par,cam]=getCameraCalibration(imloader)
file='settings/cameras.mat';
l=load(file);
val=[];
for cam=1:length(l.cameras)
    valh=imloader.gettag(l.cameras(cam).ID.tag);
    if ~isempty(val)&&strcmp(val,l.cameras(cam).ID.value)        
        val=valh;
        break
    end
end
if isempty(val)
    cam=[];
    par=[];
    return
end
partable=l.cameras(cam).par;
s=size(partable);
for k=1:s(1)
    switch partable{k,2}
        case 'fix'
            X=partable{k,3};
        case 'state dependent'
            X='0';
        case 'metadata'
            X=imloader.gettag(partable{k,4});

    end
    if ~isempty(partable{k,6})&&ischar(X)
        X=eval(partable{k,6});
    end
    par.(partable{k,1})=X;
    
end
if isempty(par.roi)
    par.roi=[0 0 par.Width par.Height];
end
end