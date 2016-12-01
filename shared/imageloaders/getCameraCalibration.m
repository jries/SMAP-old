function [par,cam]=getCameraCalibration(imloader)
file='settings/cameras.mat';
l=load(file);
for cam=1:length(l.cameras)
    val=imloader.gettag(l.cameras(cam).ID.tag);
    if ~isempty(val)
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
    whos X
    if ~isempty(partable{k,6})&&ischar(X)
        X=eval(partable{k,6});
    end
    par.(partable{k,1})=X;
    
end
end