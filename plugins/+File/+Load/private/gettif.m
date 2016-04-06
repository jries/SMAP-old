function imout=gettif(file)
imout.image=imread(file); 
sim=size(imout.image);
imout.info.Width=sim(1);
imout.info.Height=sim(2);
imout.info.roi=getRoiTif(file);
imout.info.name=file;
end