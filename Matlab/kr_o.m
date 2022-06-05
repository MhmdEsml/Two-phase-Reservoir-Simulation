function [kro, kro_s] = kr_o (sw)
% interpolation of relative permeability for water
% data set:
% swc=0.15, sor=0.3, kro_0=1.0, a=3
% kro= kro_0*(1-sw-swc)/(1-swc-sor)
kro= 1.0* (1.0-sw-0.3)^3.0/(1.0-0.15-0.3);
kro_s= 1.0*2.0* (1.0-sw-0.3)^2.0/(1.0-0.15-0.3);
if kro<0
    kro=0;
    kro_s=0;
end

