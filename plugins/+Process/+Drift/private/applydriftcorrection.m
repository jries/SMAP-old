function poso=applydriftcorrection(drift,pos)
indg=pos.frame>0;
if isfield(drift,'x') 
    poso.xnm=pos.xnm;
    poso.xnm(indg)=pos.xnm(indg)-drift.x(pos.frame(indg));
end
if isfield(drift,'y')
    poso.ynm=pos.ynm;
    poso.ynm(indg)=pos.ynm(indg)-drift.y(pos.frame(indg));
end
if isfield(drift,'z')&&isfield(pos,'znm')
    poso.znm=pos.znm;
    poso.znm(indg)=pos.znm(indg)-drift.z(pos.frame(indg));
end