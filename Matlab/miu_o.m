function Viso = miu_o (Po)
% interpolation of oil viscousity
% data set:
viso=[1.16 1.164 1.167 1.172 1.177 1.181 1.185 1.190];
po=[400 1200 2000 2800 3600 4400 5200 5600];
Viso=interp1(po,viso,Po,'linear','extrap');

