function zout=correct_3Daberrations(zcorr,zin,objectivepos)
% zcorr structure returned by calibrate3Daberrations
% zin: uncorrected z (all units: nm)
% objectivepos: position of the focal plane above the coverslip in
% nanometers
% zout: correctec z-positions
zout=zin+zcorr(ones(size(zin))*objectivepos,zin);