function [par,cam,state]=getCameraCalibration(imloader,l,silent)
if nargin<3
    silent=false;
end
replacefields={'cam_pixelsize_um','pixsize'};

par=[];
cam=[];
state=[];
argin{1}=imloader;argin{3}=silent;
if nargin<2||isempty(l)
    argin{2}=[];
    file='settings/cameras.mat';
    if ~exist(file,'file')
%         l=[];
        [par,cam,state]=askforcameramanager(imloader,'camera calibration file settings/camera.mat not found. ',silent,argin);
        return
%         [par,cam]=getCameraCalibration(imloader,[],silent);
    end
    l=load(file);
    argin{2}=l;
else
    argin{2}=l;
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
    cam=[];
    
    if ~silent
        answ=questdlg('Camera not recognized. Use default camera?');
    end
    if silent||strcmp(answ,'Yes')  
        camnames=getFieldAsVector(l.cameras,'ID','name');
        cam=find(strcmp(camnames,'Default'));
        state=1;
        if isempty(cam)&~silent
            errordlg('create Default camera with Camera Manager')
        end    
    else    
        [par,cam,state]=askforcameramanager(imloader,'Camera not recognized. ',silent,argin);
        return;
    end
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
                found(k)=false;
                break
            end
                
        end
    end
    state=find(found,1,'first');  
    if isempty(state)
%         askforcameramanager(imloader,'State of the camera could not be determined. Please use the CameraManager to define proper state. Create new state with Camera Manager now?',silent)
            [par,cam,state]=askforcameramanager(imloader,'State of the camera could not be determined. Please use the CameraManager to define proper state.',silent,argin);
            if ~isempty(par)
                return;
            end
    end
end
for k=1:s(1)
    switch partable{k,2}
        case 'fix'
            X=partable{k,3};
        case 'state dependent'
            if isempty(state)
                X=[];
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
for k=1:size(replacefields,1)
    par.(replacefields{k,2})=par.(replacefields{k,1});
end

function [paro,camo,stateo]=askforcameramanager(imloader,message,silent,argin)
paro=par;camo=cam;stateo=state;
if silent
    disp(message)
    return
end
answ=warndlg(message,'Open Camera Manager now?');
if strcmp(answ,'Yes')
    disp('please open in camera manager');
    
%     camm=CameraManager;
%     camm.imloader=imloader;
%     camm.loadimages;
%     waitfor(camm.handle)
%     argin{2}=[];
%     [paro,camo,stateo]=getCameraCalibration(argin{:});
end
end
end