classdef AnalyzeSPT<interfaces.DialogProcessor
    %
    methods
        function obj=AnalyzeSPT(varargin)        
            obj@interfaces.DialogProcessor(varargin{:}) ;
            obj.inputParameters={'sr_pixrec','numberOfLayers','sr_pos','sr_size','layers','sr_layerson'};
            obj.showresults=true;
        end
        function out=run(obj,p)
            
            out=analyzei(obj,p);
        end
        function pard=guidef(obj)
            pard=guidef;
        end
    end
end
function out=analyzei(obj,p)
[locs,indin]=obj.locData.getloc({'xnm','ynm','znm','frame','track_id'},'layer',1,'position','roi','grouping','ungrouped');
intrack=find(locs.track_id>0);
trackid=locs.track_id(intrack);
%maybe filter with minimum length here: histogram 1:N
hc=histcounts(trackid,1:max(trackid)+1);
minlen=10;
goodids=(hc>minlen);
longtrack=goodids(trackid);

loct.x=locs.xnm(intrack(longtrack));
loct.y=locs.ynm(intrack(longtrack));
loct.frame=locs.frame(intrack(longtrack));
loct.track_id=locs.track_id(intrack(longtrack));
if ~isempty(locs.znm)
    loct.z=locs.znm(intrack(longtrack));
end


ax=obj.initaxis('tracks overlay');
img=obj.getPar('sr_image');
imagesc(ax,img.rangex*1000,img.rangey*1000,img.image)
p.timediff=30;%XXX
p.mintracklength=minlen;

hold on
plottracks(loct,p,2)
hold off
ax=obj.initaxis('MSD hist');
stat=plottracks(loct,p,0);
hist((stat.slo),30)
xlabel('(slope)')

avim=sum(img.image,3);
c = max( avim(:) );
avim = avim / c;
avim2 = zeros([size(avim) 3]);
avim2(:,:,1) = avim;
avim2(:,:,2) = avim;
avim2(:,:,3) = avim;
plotTrackDiffusion(loct,p,img.rangex*1000,img.rangey*1000,avim2,-2,2);

%somewhere: write diffusion coefficient to localization data.
end
            

function out=plottracks(locs,p,show)

	pfound = max(locs.track_id);
	cols = jet(pfound); 		% default colormap
	off = zeros(pfound,1);		% double vector (all zeros)
	slo = off;
	lent = off;
	indg = false(size(locs.track_id));			% boolean vector (all false)
	ind1 = 1;
	st = length(locs.track_id);
	dt = p.timediff/1000;
	
    p.minlenplot=10;
	to=zeros(st,1);			% add a vector with zeros

	for k = 1:pfound
		
		% increment ind1 as long as
		%  - it is greater then the number of rows in tracks
		%  - the corrisponding row has track-id smaller k
		while( ind1 < st(1) && locs.track_id(ind1) < k )
			ind1 = ind1 + 1;
		end
		
		ind2 = ind1;
		
		% increment ind2 as long as
		%  - it is greater then the number of rows in tracks
		%  - the corrisponding row has track-id equal k
		while( ind2 < st(1) && locs.track_id(ind2) == k)
			ind2=ind2+1;
		end
		
		% take all rows between ind1 and ind2 (-> index always k)
		ind = ind1:ind2 - 1;
		
		% par.analv1 == min L (GUI)
		if ( length(ind) >= p.minlenplot)
			indg(k) = true; 					% mark row as checked
			x = locs.x(ind) ;
			y = locs.y(ind);
			d = zeros( length(x)-1 ,1);			% vector with length of x
			
			% track step size considered for diffusion
			mint = 2;
			maxt = 4;
			
			% d contains the sum of mean squared distances for dt = 2 : 4 (time step size)
			for l = mint:maxt
				d(l) = mean( (x(1:end-l)-x(l+1:end)).^2 + (y(1:end-l)-y(l+1:end)).^2 );
			end
			
			t = dt * (mint:maxt)';
			
			% X(:,1) = 1, X(:,2) = t
			X = [ones(size(t)) t];
			
			% X  *  coeffs  =  msd
			coeffs = X \ d(mint:maxt);
			slo(k) = coeffs(2)/1e6;         % slope  (k)
			off(k) = coeffs(1);         % offset (d)
			lent(k) = length(ind);		% track length
			
			to(ind,end) = coeffs(2);	% add slope to tracks
			
			if (show == 2)	%plot tracks
				plot(x,y,'Color',cols(k,:))
			end
			
			if (show == 1)						%plot msd vs time
				for l = 1:length(x)-1
					d(l) = mean( (x(1:end-l)-x(l+1:end)).^2 + (y(1:end-l)-y(l+1:end)).^2 );
				end
				plot(d)
				hold on
			end
		end
	end

	out.slo = slo(indg);
	out.off = off(indg);
	out.lent = lent(indg);
	out.tracks = to;

end

% plotTrackDiffusion
%
% Input
% - tracks   		 : [x y time id ...]
% - par	   	 		 : parameters (global variable)
% - imgx,imgy,avim2	 : background image parameters
% - minD,maxD        : min and max diffusion coefficient
function plotTrackDiffusion(locs,par,imgx,imgy,avim2,minD,maxD)

pfound = max(locs.track_id);

	st = length(locs.track_id);
	dt = par.timediff/1000;
    
    colors=jet(255);
	ind1 = 1;

	figure(33)
	clf
	% plot background image
	image(imgx,imgy,avim2);
	
	%Min
	uicontrol('Style', 'slider',...
		'Min',-8,'Max',2,'Value',minD,...
		'Position', [80 50 500 20],...
		'Callback', {@plotTrackDiffusionHelper,locs,par,imgx,imgy,avim2,minD,maxD,1});
	%Max	
	uicontrol('Style', 'slider',...
		'Min',-4,'Max',6,'Value',maxD,...
		'Position', [80 20 500 20],...
		'Callback', {@plotTrackDiffusionHelper,locs,par,imgx,imgy,avim2,minD,maxD,0});
	uicontrol('Style', 'pushbutton', 'String', 'Adjust Colors',...
		'Position', [600 20 100 20],...
		'Callback', 'colormapeditor');
	text = uicontrol('Style','text',...
		'Position',[600 45 100 20],...
		'String','D');

	hold on
	Min = inf;
	Max = -inf;
	usedD = [];
	for k = 1:pfound
		
		% increment ind1 as long as
		%  - it is greater then the number of rows in tracks
		%  - the corrisponding row has track-id smaller k
		while( ind1 < st(1) && locs.track_id(ind1) < k )
			ind1 = ind1 + 1;
		end
		
		ind2 = ind1 + 1;
		
		% increment ind2 as long as
		%  - it is greater then the number of rows in tracks
		%  - the corrisponding row has track-id equal k
		while( ind2 < st(1) && locs.track_id(ind2) == k)
			ind2=ind2+1;
		end
		
		% take all rows between ind1 and ind2 (-> index always k)
		ind = ind1:ind2 - 1;
		
		% par.analv1 == min L (GUI)
		if ( length(ind) >= par.mintracklength) 
			indg(k) = true; 					% mark row as checked
			x = locs.x(ind);
			y = locs.y(ind) ;
			d = zeros( length(x)-1 ,1);			% vector with length of x
			
			% track step size considered for diffusion
			mint = 2;
			maxt = 4;
			
			% d contains the sum of mean squared distances for dt = 2 : 4 (time step size)
			for (l = mint:maxt)
				d(l) = mean( (x(1:end-l)-x(l+1:end)).^2 + (y(1:end-l)-y(l+1:end)).^2 )/1e6;
			end
			
			t = dt * (mint:maxt)';
			
			% X(:,1) = 1, X(:,2) = t
			X = [ones(size(t)) t];
			
			% X  *  coeffs  =  msd
			coeffs = X \ d(mint:maxt);
            coeffs(imag(coeffs)~=0)=0;
% log10(coeffs(2) / 4)
			colormap jet;
			if coeffs(2)>0&&( log10(coeffs(2) / 4) < maxD && log10(coeffs(2) / 4) > minD )
                usedD = [usedD coeffs(2)/4];
				Min = min(coeffs(2)/4,Min);
				Max = max(coeffs(2)/4,Max);
                colcoeff=coeffs(2)/4;colcoeff(colcoeff<0)=0;
                colcoeff=(log10(colcoeff)-minD)/(maxD-minD);
                c=round(colcoeff*255); c(c<1)=1;c(c>255)=255;
                l = plot(x,y,'Color',colors(c,:));
% 				l = color_line(y,x,repmat(coeffs(2)/4,1,length(x)));
				set(l,'ButtonDownFcn',{@plotTrackD coeffs(2)/4 text})
			end
		end
    end
    
	s = (maxD - minD) / 10;
    cr=0:0.1:1;
	h = colorbar;
	set(h,'YTick',cr)
    
% 	set(h,'YLim',[minD maxD])
	set(h,'YTickLabel',{(minD:s:maxD)})
% % 	set(h,'YTickMode','auto')
% 	set(gcf, 'Renderer', 'opengl')

	hold off
	axis equal
	axis([ min(locs.x) max(locs.x) min(locs.y) max(locs.y) ])
	title(['D interval = [' num2str(minD,2) ',' num2str(maxD,2) '] - Median D = ' num2str(median(usedD),3) ' µm^2/s'])
end	
	
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% helper functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% plotTrackDiffusionHelper
% updates min and max value and executes plotTrackDiffusion
%
% Input
% - hObj	 	     : slider handle
% - event	 		 : event data
% - tracks   		 : [x y time id ...]
% - par	   	 		 : parameters (global variable)
% - imgx,imgy,avim2	 : background image parameters
% - minVal, maxVal   : time lag interval for MSD
% - newMin           : boolean (1 == new value is new min, 0 == new value is new max)
function plotTrackDiffusionHelper(hObj,event,tracks,par,imgx,imgy,avim2,minD,maxD,newMin)
	if hObj == 0
		plotTrackDiffusion(tracks,par,imgx,imgy,avim2,minD,maxD);
	else 
	    val = get(hObj,'Value');
	    if (newMin)
	        minD = val;
	    else
	        maxD = val;
	    end
	end

    plotTrackDiffusion(tracks,par,imgx,imgy,avim2,minD,maxD)
end
% plotTrackD

% Input
% - hObj : line handle
% - event: event data
% - d	 : diffusion coefficient
% - hTxt : text handle
function plotTrackD(hObj,event,d,hTxt)
	for h = findobj('LineWidth',2)
		set(h,'LineWidth',.5);
	end
	set(hObj,'LineWidth',2);
	set(hTxt,'String',strcat('D = ', num2str(d)));
end


function pard=guidef

pard.analysismode.object=struct('String',{{'statistics','interactive track exploration','diffusion maps','grid based diffusion coefficients'}},'Style','popupmenu');
pard.analysismode.position=[1,1];
pard.analysismode.Width=1.5;

pard.lentrackst.object=struct('String','Minimum length of tracks','Style','text');
pard.lentrackst.position=[2,1];
pard.lentrackst.Width=1.5;
pard.lentrackst.TooltipString=sprintf(['set this keyword to eliminate all trajectories with \n'...
            ' fewer than param.good valid positions.  This is useful \n'...
           'due to blinking noise particles in the data stream.']);
pard.lentracks.object=struct('String','2','Style','edit');
pard.lentracks.position=[2,2.5];
pard.lentracks.Width=.5;
pard.lentracks.TooltipString=pard.goodt.TooltipString;

pard.plugininfo.description=sprintf('AnalyzeSPT');
pard.plugininfo.type='ProcessorPlugin';
pard.plugininfo.Name='AnalyzeSPT';
end

