function calibrate3D(p)
% p.filelist
% p.outputfile
% p.dz
% p.modality
% p.zcorr
% p.ROIxy
% p.ROIz
% p.smoothxy
% p.smoothz
% p.gaussrange
% p.filter;
% p.zcorrframes
%p.gaussroi

%get bead positions
disp('get bead positions')
[beads,p]=images2beads_so(p);
f=figure;
p.tabgroup=uitabgroup(f);
p.ploton=false;
if contains(p.modality,'astig')
    %determine sx,sy
    disp('fit beads to get sx,sy')
    for k=1:length(beads)
        stackh=single(beads(k).stack.image);
        s=size(stackh); 
        d=round((s(1)-p.gaussroi)/2);
        stack=stackh(d+1:end-d,d+1:end-d,:);
        P=callYimingFitter(stack,1,100,4,0,0);
        beads(k).loc.PSFxpix=P(:,5);
        beads(k).loc.PSFypix=P(:,6);
        beads(k).loc.phot=P(:,3);
        beads(k).loc.bg=P(:,4);
        beads(k).f0=stackas2z_so(beads(k).loc.PSFxpix,beads(k).loc.PSFypix,beads(k).loc.frames,beads(k).loc.phot,p);
        ind=find(beads(k).loc.frames<=beads(k).f0,1,'last');
        if isnan(beads(k).f0)||isempty(ind)
            ind=1;
        end
        beads(k).psfx0=beads(k).loc.PSFxpix(ind);
        beads(k).psfy0=beads(k).loc.PSFypix(ind);
    end
    badind=isnan([beads(:).f0]);
    beads(badind)=[];
else
    f0g=round(size(beads(1).stack,3)/2);
    for k=1:length(beads)
        beads(k).f0=f0g;
    end
end
if contains(p.modality,'astig')
    [curvecal,indgoodc]=getgausscal_so(beads,p); 
else
    indgoodc=true(size(beads));
end
[stackcal,indgoods]=getstackcal_so(beads(indgoodc),p);

end

        function runold
                    
            for X=1:length(p.Xrange)-1
                for Y=1:length(p.Yrange)-1
                    for iter=1:redo
                        if p.uselocs
                            [curvecal,indgoodc]=getcurvecal(beadh,p,X,Y,axall,indgood); 
                        else
                            indgoodc=indgood;
                        end
                        if  p.fitcsplinec
                             disp('get stack spline calibration')
                            [stackcal,indgoods]=getstackcal(beadh,p,X,Y,axall,indgood);
                            if p.uselocs
                                for Z=length(curvecal):-1:1
                                    allcal(Z)=copyfields(curvecal(Z),stackcal(Z),{'splinefit'});
                                end
                            else
                                allcal=stackcal;
                            end
                        else
                            allcal=curvecal;
                            indgoods=indgoodc;
                        end

                        for Z=1:length(indgoods)
                            indgood{Z}=indgoods{Z}&indgoodc{Z};
                        end
    %                     indhf=find(indh);
    %                     indgood(indhf)=indgood(indhf)&indgoods&indgoodc;
    %                     indgood=indgood&indgoods&indgoodc;
                        SXY(X,Y,:)=allcal;
                    end
                end       
            end
            plotcurves(obj,SXY,axall,p)
            %save
            lastf=obj.getPar('lastSMLFile');
            if ~isempty(lastf)&&p.beadsource.Value~=3
            [path,file,ext]=fileparts(lastf);
            file=[file ext];
            path=[path filesep];
            else
                [file,path]=uiputfile([pathhere '3d_3Dcal.mat']);
                if ~file
                    error('no file selected')
                end
            end
            if p.spatialcalibration && p.zcalc
                adds=p.zfilter.selection;
            else
                adds='';
            end
            file=strrep(file,'_sml',['_3Dcal' adds]);
            save([(path)  file],'SXY')
        end



function axall=getaxes(p)
            ax=initaxis(p.resultstabgroup,'beads');
            axall.htbeads=uitabgroup(ax.Parent);
            axall.axbeads=maketgax(axall.htbeads,'scatter');   
            
            if p.fitcsplinec
                ax=initaxis(p.resultstabgroup,'spline fit');
                axall.allsplines=uitabgroup(ax.Parent);
                ax=maketgax(axall.allsplines,'scatter');
%                 ax=initaxis(axall.allsplines,'scatter');
                axall.hspline_scatter=uitabgroup(ax.Parent);
                
                ax=maketgax(axall.allsplines,'PSF');
                axall.hspline_psf=uitabgroup(ax.Parent);
                
                ax=maketgax(axall.allsplines,'Overlay');
                axall.hspline_overlay=uitabgroup(ax.Parent);
                
                 ax=maketgax(axall.allsplines,'PSFz');
                axall.hspline_psfz=uitabgroup(ax.Parent);
                
                ax=maketgax(axall.allsplines,'PSFx');
                axall.hspline_psfx=uitabgroup(ax.Parent);
                
                ax=maketgax(axall.allsplines,'validate');
                axall.hspline_validate=uitabgroup(ax.Parent);
                
                ax=maketgax(axall.allsplines,'stripes');
                axall.hspline_stripes=uitabgroup(ax.Parent);
            end  
                
            %init validation and summary axes
            ax=initaxis(p.resultstabgroup,'splines');
            axall.htsplines=uitabgroup(ax.Parent);
            axall.axsplines=maketgax(axall.htsplines,'summary'); 

            if p.fitzsxsyc
                ax=initaxis(p.resultstabgroup,'sx^2-sy^2');
                axall.htsx2sy2=uitabgroup(ax.Parent);
                axall.axsxsys=maketgax(axall.htsx2sy2,'summary');  
                axall.axsxs22=maketgax(axall.htsx2sy2,'validation'); 
            end
            if p.fitzc
                ax=initaxis(p.resultstabgroup,'z fit');
                axall.htzfit=uitabgroup(ax.Parent);
                axall.axzfits=maketgax(axall.htzfit,'summary');
            end
            if p.calculateZSxSy
                ax=initaxis(p.resultstabgroup,'Z(Sx,Sy)');
                axall.hzsx=uitabgroup(ax.Parent);
 
                axall.axzlut=maketgax(axall.hzsx,'validation'); 
            end
              
end



       function plotcurves(obj,SXY,axall,p)
            %plot results
            axes(axall.axbeads)
            sp=SXY(:);
            legends={SXY(:).legend};
            for k=1:numel(sp)
                if ~isempty(sp(k).curve)
                    xpos=vertcat(sp(k).curve(:).xpos);
                    ypos=vertcat(sp(k).curve(:).ypos);
                    plot(xpos,ypos,'k.');hold on;
                    text(mean(sp(k).Xrange),mean(sp(k).Yrange),sp(k).legend)
                end
            end
            
            for k=1:length(p.Xrange)
                line([p.Xrange(k),p.Xrange(k)],[p.Yrange(1) p.Yrange(end)])
            end
            for k=1:length(p.Yrange)
                line([p.Xrange(1),p.Xrange(end)],[p.Yrange(k) p.Yrange(k)])
            end
            
            
            sp=[SXY(:).spline];
%             zt=zrangeall(1):0.01:zrangeall(2);
%             z0=zf1-zpos;
%             z0=zshift;
            if ~isempty(sp)
                z0=0;
                zt=z0+p.zrangeuse(1):0.01:z0+p.zrangeuse(2);

                linecolors=lines(length(sp));
                axes(axall.axsplines)
                hold off;

                for k=1:numel(sp)
                    pl2(k)=plot(axall.axsplines,zt,sp(k).x(zt),'Color',linecolors(k,:));
                    hold on;
                    plot(axall.axsplines,zt,sp(k).y(zt),'Color',linecolors(k,:));

                end
                ylim(axall.axsplines,[0 5])
                legend(pl2,legends);
            end
            
            if ~isempty(SXY(1).Sx2_Sy2)
                sp={SXY(:).Sx2_Sy2};
                s=-30:0.5:30;
                axes(axall.axsxsys)
                hold off;
                for k=1:numel(sp)
                    if ~isempty(sp{k})
                    plot(axall.axsxsys,sp{k}(s),s);hold on;
                    end
                end
                xlabel(axall.axsxsys,'z')
                ylabel(axall.axsxsys,'Sx^2-Sy^2')
                xlim(axall.axsxsys,p.zrangeuse)
                legend(legends);
            end
            
            if ~isempty(SXY(1).fitzpar)
                sp={SXY(:).fitzpar};


                axes(axall.axzfits)
                hold off
                xlabel(axall.axzfits,'z')
                ylabel(axall.axzfits,'Sx, Sy')
                ylim(axall.axzfits,[0 5])

                hold off;
    %             nleg={};
                pl=[];
                for k=1:numel(sp)
                    if ~isempty(sp{k})
                        [sx,sy]=getsxfromzfitpar(zt,sp{k},z0); 
                        pl(k)=plot(axall.axzfits,zt,sx,'Color',linecolors(k,:));
                         hold on
                        plot(axall.axzfits,zt,sy,'Color',linecolors(k,:))

    %                     nleg{k}=num2str(k);
                    end
                end

                legend(pl,legends);
            end
            %cross check and validate with bead positions
%             axes(axzlut);
            
%             hold off
            zt=p.zrangeuse(1):p.dz:p.zrangeuse(2);
             
            for k=1:numel(SXY)
                if ~isempty(SXY(k).curve)
                    sxa=vertcat(SXY(k).curve(:).sx);
                    sya=vertcat(SXY(k).curve(:).sy);
                    za=vertcat(SXY(k).curve(:).z);

                    if ~isempty(SXY(k).splineLUT)
                        zha=zfromSXSYLut(SXY(k).splineLUT,sxa,sya);
                        dzba=bindata(za,za-zha,zt,'mean');
                        dzs=bindata(za,za-zha,zt,'std');
                        plot(axall.axzlut,za,za-zha,'.','MarkerSize',2)
                        axall.axzlut.NextPlot='add';
                        plot(axall.axzlut,zt,dzba,'k',zt,dzba+dzs,'k',zt,dzba-dzs,'k')
                        plot(axall.axzlut,p.zrangeuse,[0 0],'k');
                         ylim(axall.axzlut,yrange)
                         xlim(axall.axzlut,p.zrangeuse)
                    end
                    if ~isempty(SXY(k).Sx2_Sy2)
                        zha2=zfromSx2_Ss2(SXY(k).Sx2_Sy2,sxa,sya);
                         plot(axall.axsxs22,za,za-zha2,'.','MarkerSize',2)
                         axall.axsxs22.NextPlot='add';
                        dzba2=bindata(za,za-zha2,zt,'mean');
                        dzs2=bindata(za,za-zha2,zt,'std');
                        plot(axall.axsxs22,zt,dzba2,'k',zt,dzba2+dzs2,'k',zt,dzba2-dzs2,'k')
                    end
                end
             end
            
            if isfield(axall,'axsxs22')
                plot(axall.axsxs22,p.zrangeuse,[0 0],'k');
                yrange=[-100 100];

                ylim(axall.axsxs22,yrange)
                xlim(axall.axsxs22,p.zrangeuse)   
            end
            obj.SXY=SXY;  
       end


function splineLUT=cal_splineLUT(spline,p)
srange=p.Smin:p.Sd:p.Smax;
splineLUT=gets2z(spline,srange);
% splineLUT=[];
end


function getcoords(a,b,obj)
locsall=obj.locData.getloc({'frame','xnm','ynm','PSFxnm','PSFynm','filenumber','phot'},'position','all','layer',1,'removeFilter','filenumber');
p=obj.getAllParameters;
x=locsall.xnm/p.cam_pixelsize_nm;
y=locsall.ynm/p.cam_pixelsize_nm;
img=myhist2(x,y,1,1,[0 512.5],[0 512.5]);
f1=figure;
imagesc(img');
h=imrect;
pos=wait(h);% [x y wx wy]
delete(f1);
answ=inputdlg({'number of rows','number of columns'},'set tiles', 1,{'2','1'});
if isempty(answ)
    return;
end
nx=str2double(answ{2});
ny=str2double(answ{1});

% pos=pos([2 1 4 3]);
p.Xmin=round(pos(1)); p.Ymin=round(pos(2)); p.Xmax=round(pos(1)+pos(3)); p.Ymax=round(pos(2)+pos(4));
p.Xd=floor(pos(3)/nx);
p.Yd=floor(pos(4)/ny);
obj.setGuiParameters(p)

end

