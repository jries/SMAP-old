function lo=saveaspix(locData,varargin)
if any(contains(varargin,'grouped'))
    ll=locData.grouploc;
else
    ll=locData.loc;
end

if any(contains(varargin,'nm'))
    cp=[1 1];
    roi=[0 0];
else
cp=locData.files.file(1).info.cam_pixelsize_um*1000;
roi=locData.files.file(1).info.roi;
end
try
    RIM=locData.history{1}.children.fitparamters.fitterGUI.children.MLE_GPU_Yiming.refractive_index_mismatch;
catch err
    disp('could not find history')
    RIM=0.8;
end
lo=copyfields([],ll,{'phot','bg','filenumber','frame','numberInGroup'});
lo.x=ll.xnm/cp(1)-roi(1);
lo.y=ll.ynm/cp(2)-roi(2);
lo.z=ll.znm/RIM;
lo.locprec=ll.locprecnm/cp(1);
lo.locprecz=ll.locprecznm/RIM;
lt=struct2table(lo);
[f,p]=uiputfile('*.csv');
writetable(lt,[p f]);
