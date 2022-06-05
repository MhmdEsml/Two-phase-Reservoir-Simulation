function [Phi,Phi_p] =poro(Po,phi0)

Phi = phi0 * exp(5 * 10^-6 * (Po-5000));
Phi_p = phi0 * 5 * 10^-6 * exp(5 * 10^-6 * (Po-5000));
