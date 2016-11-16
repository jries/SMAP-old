classdef TifSaver<interfaces.DialogProcessor
    methods
        function obj=TifSaver(varargin)  
                obj@interfaces.DialogProcessor(varargin{:}) ;
                obj.inputParameters={'filelist_long','mainfile','mainGui','numberOfLayers','sr_layerson','sr_pixrec','layers','sr_image','sr_pos'};
        end
        
        function out=save(obj,p)
            p2=obj.getGuiParameters;
            obj.status('save tif file')
            fn=p.filelist_long.selection;
            [path,file]=fileparts(fn);
            of=[path filesep file p2.img_ext.selection];
            ind=2;
            while exist(of,'file')
                of=[path filesep file '_' num2str(ind) p2.img_ext.selection];
                ind=ind+1;
            end
            txt=filterpar(obj,p);
            serchstr={['*' p2.img_ext.selection];['*' strjoin(p2.img_ext.String,';*')]};
            [f,path]=uiputfile(serchstr,'select output file for image', of);
            if f
                img=obj.getPar('sr_image');
                res=ones(2,1)/p.sr_pixrec*2.5e4;       
                description=sprintf(txt);
                [~,~,ext]=fileparts(f);
                switch ext
                    case '.tif'
                        imwrite(img.image,[path f],'Description',description,'Resolution',res);
                    case '.png'
                        imwrite(img.image,[path f],'Description',description,'XResolution',res(1),'YResolution',res(2));
                end
            end
            obj.status('save done')
          
        end
        function pard=guidef(obj)
           pard.plugininfo.type='SaverPlugin';
           
            pard.img_ext.object=struct('Style','popupmenu','String',{{'.tif','.png'}});
            pard.img_ext.position=[1,1];
            pard.img_ext.Width=2;
        end
        function run(obj,p)
            obj.save(p)
        end        

    end
end

function txt=filterpar(obj,p)
txt='SMAP \n';
txt=[txt 'pixelsize(nm) \t' num2str(p.sr_pixrec) '\n'];
txt=[txt 'position (nm) \t' num2str(p.sr_pos (1:2)) '\n'];
for k=1:length(p.sr_layerson)
    if p.sr_layerson(k)
        lp=['layer' num2str(k)];
        txt=[txt lp ':\n'];
        txt=[txt p.([lp '_']).ch_filelist.selection '\n'];
        
        filn=p.([lp '_']).ch_filelist.Value;
        txt=[txt p.filelist_long.String{filn} '\n'];
        try
            fn=obj.locData.files.file(filn).info.filename;
            fn=strrep(fn,'\','/');
            txt=[txt fn '\n'];
        catch
        end
        txt=[txt 'channels: ' num2str(p.([lp '_']).channels) '\n'];
        txt=[txt 'grouping: ' num2str(p.([lp '_']).groupcheck) '\n'];
        txt=[txt 'quantile/Imax: ' num2str(p.([lp '_']).imax_min) '\n'];
        txt=[txt 'color range: \t' num2str(p.([lp '_']).colorfield_min) ' : \t' num2str(p.([lp '_']).colorfield_max) '\n'];
        txt=[txt 'remove outside c-range: ' num2str(p.([lp '_']).remout) '\n'];
        if strcmp(p.([lp '_']).renderfield.selection,'field')
            txt=[txt 'render field: ' (p.([lp '_']).render_colormode.selection) '\n'];
        end
        if strcmp(p.([lp '_']).renderfield.selection,'z')
            txt=[txt 'render field: znm'  '\n'];
        end
        s=obj.getPar(['layer' num2str(k) '_filtertable']);
        ss=size(s);
        txt=[txt 'filters:'  '\n'];
        for l=1:ss(1)
            if s{l,7} &&~strcmp(s{l,1},'xnm')&&~strcmp(s{l,1},'ynm')&&~strcmp(s{l,1},'filenumber')
                line=s(l,1:6);
                th=[line{1} ':\t' num2str(line{2}) ' : \t' num2str(line{6}) '\n'];
                txt=[txt th];
            end
            
        end
%         disp(sprintf(txt));
    end
    
end
end