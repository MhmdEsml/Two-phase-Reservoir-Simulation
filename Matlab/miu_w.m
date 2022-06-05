function Visw = miu_w (Po)
% interpolation of water viscousity
% data set:
visw=[0.52341 0.56341];
po=[3600 3900];
Visw=interp1(po,visw,Po,'linear','extrap');


