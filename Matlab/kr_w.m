function [krw, krw_s] = kr_w (sw)
% interpolation of relative permeability for oil
% data set:

% swc=0.15, sor=0.3, krw_0=0.5, a=3
% krw= krw_0*(1-sw-swc)/(1-swc-sor)
krw= 0.5* (sw-0.15)^3.0/(1.0-0.15-0.3);
krw_s= -0.5*3* (sw-0.15)^2.0/(1.0-0.15-0.3);
if krw<0
    krw=0;
    krw_s=0;
end
