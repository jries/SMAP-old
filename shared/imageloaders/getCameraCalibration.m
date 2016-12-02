function [par,cam]=getCameraCalibration(imloader,l,silent)
if nargin<3
    silent=false;
end
if nargin<2||isempty(l)
    file='settings/cameras.mat';
    if ~exist(file,'file')
        askforcameramanager(imloader,'camera calibration file settings/camera.mat not found. Create new file with Camera Manager?',silent)
        par=[];
        cam=[];
    end
    l=load(file);
end
val=[];
for cam=1:length(l.cameras)
    valh=imloader.gettag(l.cameras(cam).ID.tag);
    if ~isempty(valh)&&strcmp(valh,l.cameras(cam).ID.value)        
        val=valh;
        break
    end
end
if isempty(val)
    askforcameramanager(imloader,'Camera not recognized. Create new camera with Camera Manager?',silent)
    cam=[];
    par=[];
    return
end
partable=l.cameras(cam).par;
s=size(partable);
if ~isempty(strcmp(partable(:,2),'state dependent'))
    %get state
%     found=false;
    for k=length(l.cameras(cam).state):-1:1
        deftab=l.cameras(cam).state(k).defpar;
        found(k)=true;
        for k2=1:size(deftab,1)
            if strcmp('select',deftab(k2,1) )|| isempty(deftab{k2,2})
                continue
            end
            valh=imloader.gettag(deftab{k2,1});
            if ~strcmp(valh,deftab{k2,2})
                found(l)=false;
                break
            end
                
        end
    end
    state=find(found);  
    if isempty(state)
        askforcameramanager(imloader,'State of the camera could not be determined. Please use the CameraManager to define proper state. Create new state with Camera Manager now?',silent)
    end
end
for k=1:s(1)
    switch partable{k,2}
        case 'fix'
            X=partable{k,3};
        case 'state dependent'
            if isempty(state)
                X='';
            else
                X=l.cameras(cam).state(state).par{k,2}; 
            end
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

function askforcameramanager(imloader,message,silent)
if silent
    disp(message)
    return
end
answ=questdlg(message,'Open Camera Manager now?');
if strcmp(answ,'Yes')
    camm=CameraManager;
    camm.imloader=imloader;
    camm.loadimages;
end
end