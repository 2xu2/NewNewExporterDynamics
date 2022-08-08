function [Pi] = profit(X, xi, Evalue, Qvalue, theta, Cstar, alphan, alphak, w, r, tau, P, Pstar)
% Profit function for a plant in new exporter dynamics of Ruhl, Willis (2013)
% X is the export decision. X = 0 and 1
% Qvalue and Evalue are realizations of simulated AR(1) processes (vectors)
% theta is the elasticity of substitution
% Cstar is the size of world aggregate demand relative to domestic demand
% alphan and alphak are parameters of the Copp-Douglas production functio
% w is wage, r is the rent price of capital 
% State space is X, Epsilon and Q
% xi is the iceberg cost. Note when X = 0, value of xi doesn't matter
% P and Pstar are aggragate price levels in domestic and foreign markets
% tau is the ad-valorem tariff rate

% Not fastest, but clear
M = (1 + X * xi^(1 - theta) * Cstar  * (1 - tau)^theta * ((Qvalue.*(Pstar/P)).^theta)).^(1/theta) .* P .* Evalue;
%define two cobb douglas coefficients
a = alphan * (theta - 1)/theta;
b = alphak * (theta - 1)/theta;
Pi = (1-a-b) * (a/w)^(a/(1-a-b)) * (b/r)^(b/(1-a-b)) .* M.^(1/(1-a-b));
end


