classdef NPCgeomtryQuantify<interfaces.SEEvaluationProcessor
    properties
        savedevals
    end
    methods
        function obj=NPCgeomtryQuantify(varargin)        
                obj@interfaces.SEEvaluationProcessor(varargin{:});
        end
        function makeGui(obj)
            makeGui@interfaces.SEEvaluationProcessor(obj);
%             obj.guihandles.saveimagesb.Callback={@saveimagesb_callback,obj};
        end
        function out=run(obj,p)
            try
            out=runintern(obj,p);
            catch err
                err
                out=[];
            end
         
        end
        
        function pard=guidef(obj)
            pard=guidef;
        end
    end
end

function pard=guidef
pard.Rt.object=struct('Style','text','String','R (nm):');
pard.Rt.position=[1,1];
pard.Rt.Width=1;

pard.R.object=struct('Style','edit','String','50');
pard.R.position=[1,2];
pard.R.Width=1;

pard.dRt.object=struct('Style','text','String','dR:');
pard.dRt.position=[1,3];
pard.dRt.Width=1;

pard.dR.object=struct('Style','edit','String','20');
pard.dR.position=[1,4];
pard.dR.Width=1;

pard.plugininfo.type='ROI_Evaluate';
pard.inputParameters={'numberOfLayers','sr_layerson','se_cellfov','se_sitefov','se_siteroi','layer1_','layer2_','se_sitepixelsize'};
end


function out=runintern(obj,p)
locs=obj.getLocs({'xnm','ynm','znm','locprecznm','locprecnm'},'layer',1,'size',p.se_siteroi);
dz=5;
z=-200:dz:200;
hz=hist(locs.znm-obj.site.pos(3),z);
% hz=hz-mean(hz);
ac=myxcorr(hz,hz);
% if obj.display
ax1=obj.setoutput('profile');

fitresult=createFit(z, hz,ax1)
title(ax1,fitresult.d)
% plot(ax1,z,hz)
ax2=obj.setoutput('correlation');
plot(ax2,z,ac)


% end
out.ac=ac;
out.dz=dz;
out.hist=hz;
out.z=z;
out.Gaussfit=fitresult;
end


function [fitresult, gof] = createFit(z, hz,ax)

[xData, yData] = prepareCurveData( z, hz );

% Set up fittype and options.
ft = fittype( 'a1*exp(-((x-b)/c)^2) + a2*exp(-((x-b-d)/c)^2)', 'independent', 'x', 'dependent', 'y' );
opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
opts.Display = 'Off';
opts.Lower = [0 0 -Inf 0 -Inf];
opts.StartPoint = [max(hz) max(hz) -50 20 100 ];

% Fit model to data.
[fitresult, gof] = fit( xData, yData, ft, opts );

axes(ax)
h = plot(fitresult, xData, yData);


legend( h, 'hz vs. z', 'untitled fit 1', 'Location', 'NorthEast' );
% Label axes
xlabel z
ylabel hz
grid on

end
