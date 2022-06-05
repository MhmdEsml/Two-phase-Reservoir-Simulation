function [bo, inv_bo_p]=b_o(p)
% oil formation factor:
bo = 0.9802 * exp( -10.0^(-5.0) * (p-3600.0));
inv_bo_p = (1/0.9802) * 10.0^(-5.0) * exp( 10.0^(-5.0) * (p-3600.0));