global se
sites=se.sites;
ls=length(sites);
s1=sites(1);
im1=s1.image.image;
sim=size(im1);
imall=zeros(sim(1),sim(2),ls,'uint8');
nlocs=zeros(ls,1);
corners=zeros(ls,1);
for k=1:ls
    sh=sites(k);
    imall(:,:,k)=uint8(sum(sh.image.image,3)*255/3);
    nlocs(k)=sh.evaluation.generalStatistics.Nlayers;
    corners(k)=sh.evaluation.NPCLabelingQuantify.numcorners;
    
end
tout=table((1:ls)',nlocs,corners);
[f,p]=uiputfile('NPCsimul.tif');
saveastiff(imall,[p f]);
% imwrite(imall,[p f]);
[~,ff]=fileparts(f);
writetable(tout,[p filesep ff '.csv']);