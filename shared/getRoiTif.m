function roi=getRoiTif(file)
p=fileparts(file);
pm=[p filesep 'metadata.txt'];
if ~exist(pm,'file')
    mf=dir([p filesep '*metadata.txt']);
    pm=[p filesep mf(1).name];
end
if exist(pm,'file')
    minfo=fileread(pm);
    searchstr='ROI":';
    ind=strfindfast(minfo,searchstr,1);
    txt1=getval(minfo,ind+length(searchstr));%,'[',']');
    txt1=strrep(txt1,'-',' ');
   
   roi=str2num(txt1); 
   if isempty(roi)
        txt2=getval(minfo,ind+length(searchstr),'[',']');
         roi=str2num(txt2);
   end

else
    roi=[0 0 512 512];
end
roi
    
    
function txt=getval(minfo,index,sep1,sep2)
if nargin==2
    sep1='"';
    sep2='"';
end

while minfo(index)~=sep1
    index=index+1;
end
i1=index+1;
index=index+1;
while minfo(index)~=sep2
    index=index+1;
end
txt=minfo(i1:index-1);