function poso=applydriftcorrection(drift,pos)
if isfield(drift,'x')
poso.xnm=pos.xnm-drift.x(pos.frame);
end
if isfield(drift,'y')
    poso.ynm=pos.ynm-drift.y(pos.frame);
end
if isfield(drift,'z')&&isfield(pos,'znm')
    poso.znm=pos.znm-drift.z(pos.frame);
end