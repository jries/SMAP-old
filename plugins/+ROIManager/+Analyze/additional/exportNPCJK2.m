global se
sites=se.sites;
ls=length(sites);
s1=sites(1);
im1=s1.image.image;
sim=size(im1);

nlocs=zeros(ls,1);
corners=zeros(ls,1);
filenumber=zeros(ls,1);
% [f,p]=uiputfile('NPCsimul.tif');
    indtab=1;
for ff=1:length(se.files)
    indim=1;
    imall=zeros(sim(1),sim(2),ls,'uint8');
    for k=1:ls
        sh=sites(k);
        if sh.info.filenumber==ff
            imall(:,:,indim)=uint8(sum(sh.image.image,3)*255/3);
            nlocs(indtab)=sh.evaluation.generalStatistics.Nlayers;
            corners(indtab)=sh.evaluation.NPCLabelingQuantify.numcorners;
            
            assignedcorners(indtab)=sh.evaluation.NPCLabelingQuantify.numbercornerassined;
            numfoundint(indtab)=sh.evaluation.NPCLabelingQuantify.numfoundint;
            cornersfiltered(indtab)=sh.evaluation.NPCLabelingQuantify.numcornersfiltered;
            numfoundrat(indtab)=sh.evaluation.NPCLabelingQuantify.numfoundrat;
            plabel(indtab)=se.files(ff).info.simulationParameters
            plabel(indtab)=
            plabel(indtab)=
            
            filenumber(indtab)=ff;
            indtab=indtab+1;
            indim=indim+1;
        end

    end
    imall(:,:,indim:end)=[];
%     saveastiff(imall,[p f]);
end


tout=table((1:ls)',nlocs,corners,filenumber);

% imwrite(imall,[p f]);
% [~,ff]=fileparts(f);
% writetable(tout,[p filesep ff '.csv']);
% 
