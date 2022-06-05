function [bw, inv_bw_p]=b_w(p)
% oil formation factor:
bw =1.029 * exp( -10^-6.0 * (p-4014.7));
inv_bw_p =  (1/1.029) * 10.0^(-6.0) * exp( 10.0^(-6.0) * (p-4014.7));