function [b,p]=images2beads_so(p)
addpath('bfmatlab')
fs=p.filtersize;
h=fspecial('gaussian',round(fs*3),fs);
fmax=0;
roisize=p.ROIxy;
roisizeh=round(1.5*(p.ROIxy-1)/2); %create extra space if we need to shift;
rsr=-roisizeh:roisizeh;
filelist=p.filelist;
b=[];
for k=1:length(filelist)
    imstack=getfile(filelist{k});
    imstack=imstack-min(imstack(:)); %fast fix for offset;
    mim=mean(imstack,3);
    mim=filter2(h,mim);
    maxima=maximumfindcall(mim);
    figure(88);
    imagesc(mim);
    drawnow
    int=maxima(:,3);
    
    mimc=mim(roisize:end-roisize,roisize:end-roisize);
    mmed=myquantile(mimc(:),0.3);
    imt=mim(mim<mmed);
    cutoff=mean(imt(:))+max(3*std(imt(:)),(max(int)-mean(imt(:)))/5);
    maxima=maxima(int>cutoff,:);
    hold on
    plot(maxima(:,1),maxima(:,2),'mo')
    hold off
    numframes=size(imstack,3);
    bind=length(b)+size(maxima,1);
%     bold=size(maxima,1);
    for l=1:size(maxima,1)
        b(bind).loc.frames=(1:numframes)';
        b(bind).loc.filenumber=zeros(numframes,1)+k;
        b(bind).filenumber=k;
        b(bind).pos=maxima(l,1:2);
        try
            b(bind).stack.image=imstack(b(bind).pos(2)+rsr,b(bind).pos(1)+rsr,:);
            b(bind).stack.framerange=1:numframes;
            b(bind).isstack=true;
            
        catch err
            b(bind).isstack=false;
            err
        end
        bind=bind-1;
    end
    fmax=max(fmax,numframes);
%     b(bold:bind)=[];
end
b=b([b(:).isstack]);

p.fminmax=[1 fmax];
p.cam_pixelsize_um=[1 1]/1000;
p.pathhere=fileparts(filelist{1});
end


function img=getfile(file)
try
imgall=bfopen(file);
catch
    disp(['could not open file ' file])
    img=[];
    return
end
sim=size(imgall{1}{1});
numframes=length(imgall{1});
img=zeros(sim(1),sim(2),numframes);
for k=1:numframes
    img(:,:,k)=imgall{1}{k};
end
end