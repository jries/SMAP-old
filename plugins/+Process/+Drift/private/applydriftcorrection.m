function poso=applydriftcorrection(drift,pos)
poso.xnm=pos.xnm-drift.x(pos.frame);
poso.ynm=pos.ynm-drift.y(pos.frame);