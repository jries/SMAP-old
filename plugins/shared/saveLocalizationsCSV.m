function  [fdesc]=saveLocalizationsCSV(locData,file,saveroi,numberOfLayers,sr_layerson)
if nargin<3
    saveroi=false;
end
if saveroi
    [~,indg]=locData.getloc('xnm','position','roi');  
    for k=1:numberOfLayers
        if sr_layerson(k)
           [~,indgh]=locData.getloc('xnm','layer',k); 
           if length(indgh)~=length(indg)
               disp('save visible works only for ungrouped data')
           end
           indg=indg&indgh;
       
        end
    end
    
else
    indg=[];
end
loc=locData.savelocs([],indg).loc; 
numlocs=length(loc.frame);
dato=zeros(numlocs,6);
dato(:,1)=1:numlocs;
dato(:,2)=loc.frame;
dato(:,3)=loc.xnm;
dato(:,4)=loc.ynm;
if isfield(loc,'znm')
    dato(:,5)=loc.znm;
end
dato(:,6)=loc.phot;

csvwrite(file,dato);
path=fileparts(file);
if 1
 %save description file
fdesc=[path filesep 'file-description.xml'];
fid=fopen(fdesc,'w');
fprintf(fid,'<?xml version="1.0" encoding="UTF-8"?> \n <description> \n');
fprintf(fid,' <firstrow>0</firstrow>\n');
fprintf(fid,' <frame>1</frame>\n');
fprintf(fid,' <xnano>2</xnano>\n');
fprintf(fid,' <ynano>3</ynano>\n');
fprintf(fid,' <znano>4</znano>\n');
fprintf(fid,' <intensity>5</intensity>\n');
fprintf(fid,' <separator>COMMA</separator>\n');
fprintf(fid,' <xshift>0</xshift>\n');
fprintf(fid,' <yshift>0</yshift>\n');
fprintf(fid,' <zshift>0</zshift>\n');
fprintf(fid,'</description>');
fclose(fid);
end
end