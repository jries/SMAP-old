classdef Phase2z4Pi<interfaces.DialogProcessor
    % sideview Side view from ROI
    methods
        function obj=Phase2z4Pi(varargin)        
            obj@interfaces.DialogProcessor(varargin{:}) ;
             obj.showresults=true;

        end
        
        function out=run(obj,p)
            fitterpath=[fileparts(obj.getPar('maindirectory')) filesep 'ries-private' filesep 'PSF4Pi'];
            addpath(fitterpath)
            locsall=obj.locData.getloc({'znm','phase','znmerr','phaseerr','frame'});
            locs=obj.locData.getloc({'znm','phase','znmerr','phaseerr','frame','filenumber'},'layer',find(obj.getPar('sr_layerson')),'position','fov');
            zastig=locs.znm;
            phase=mod(locs.phase,2*pi);
            zastigerr=locs.znmerr;
            phaseerr=locs.phaseerr;
            cal3D=obj.locData.files.file(locs.filenumber(1)).savefit.cal3D;
            frequency=cal3D.frequency/cal3D.dz;
            
            numwindows=50; windowsize=ceil(max(locs.frame)/numwindows);
            framepos=0:windowsize:max(locs.frame);
            
            z0=0;
            for k=1:length(framepos)-1
                inframe=locs.frame>framepos(k)&locs.frame<framepos(k+1);
                z0=getz0phase(zastig(inframe),phase(inframe),frequency,z0);
                z0all(k)=z0;
            end
            frameposc=framepos(1:end-1)+(framepos(2)-framepos(1))/2;
            z0int=fit(frameposc',z0all','smoothingspline');
            figure(88);plot(frameposc,z0all,frameposc,z0int(frameposc))
            z0=z0int(locsall.frame);
            zph=z_from_phi_JR(locsall.znm,mod(locsall.phase,2*pi),frequency,z0);
            obj.locData.setloc('zphase',zph);
            obj.locData.setloc('zastig',locsall.znm);
            
            % znm average:
%             zpherr=phaseerr/2/frequency;
%             obj.locData.setloc('zphasecorr',zph-z0);
            
            obj.locData.setloc('znm',zph);
            out=0;
        end
        function pard=guidef(obj)
            pard=guidef;
        end
    end
end

function pard=guidef
% pard.text1.object=struct('String','parameters','Style','text');
% pard.text1.position=[1,1];
% 
% pard.text2.object=struct('String','zmin','Style','text');
% pard.text2.position=[2,1];
% pard.text3.object=struct('String','zmax','Style','text');
% pard.text3.position=[3,1];
% 
% pard.zmin.object=struct('Style','edit','String',-400); 
% pard.zmin.position=[2,2];
% pard.zmax.object=struct('Style','edit','String',400); 
% pard.zmax.position=[3,2];
% 
% pard.pixauto.object=struct('Style','checkbox','String','set pixelsize (x,z)','Value',0);
% pard.pixauto.position=[4,1];
% pard.pixrecset.object=struct('Style','edit','String','5, 5'); 
% pard.pixrecset.position=[4,2];

% pard.plugininfo.description= 'Side view from ROI';
pard.plugininfo.type='ProcessorPlugin';
end

function z0=getz0phase(zastig,phase,frequency,z0)
% phasez=mod((zastig-z0)*2*frequency,2*pi);
% cyclicaverage(mod(phase-phasez,2*pi),2*pi)
% z0=0
zfp=phase/2/frequency;
dz=zfp-zastig+pi/frequency/2+z0;
dzm=mod(dz,pi/frequency);
z0=-cyclicaverage(dzm,pi/frequency)+pi/frequency/2+z0;

% figure(88);histogram(dzm)
% waitforbuttonpress
% phasez=mod((zastig-z0)*2*frequency,2*pi);
% dphase=phase-phasez;
% if sum(dphase>0)>sum(dphase<0)
%     dphasem=mean(dphase(dphase>0));
%     phx=-pi/frequency;
% else
%     dphasem=mean(dphase(dphase<0));
%     phx=-pi/frequency;
% end
% dz=dphasem/2/frequency;
% z0=z0+dz-phx;
% 
% phasez1b=mod((zastig-z0)*2*frequency,2*pi);

% iter=10;
% err=1;
% for k=1:iter
% phasez2=mod((zastig-z0)*2*frequency,2*pi);
% dphasez2=phase-phasez2;
% dz=mean(dphasez2)/2/frequency;
% z0=z0+dz;
% if abs(dz)<err
%     break
% end%roubst mean later?
% 
% end
% 
% phasez2=mod((zastig-z0)*2*frequency,2*pi);
% figure(88);plot(zastig,phase,'.',zastig,phasez2,'+')
% % % figure(88);plot(zastig,phase,'.',zastig,phasez,'+',zastig,phasez2,'*',zastig,phasez1b,'x')
% waitforbuttonpress
% f=@(z0,x) mod((x-z0)*2*frequency,2*pi);
% fp=fit(zastig,phase,f,'StartPoint',z0);
% phasez2=mod((zastig-fp.z0)*2*frequency,2*pi);
% figure(88);plot(zastig,phase,'.',zastig,phasez2,'+')
% waitforbuttonpress
% z0=fp.z0;
end
