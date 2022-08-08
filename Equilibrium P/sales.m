function [domestic,export] = sales(X, xi, Evalue, Qvalue, theta, Cstar, alphan, alphak, w, r, tau, P, Pstar)
% Domestic and export sales for firm decision

% Not fastest, but clear
M = (1 + X * xi^(1 - theta) * (1 - tau)^theta * Cstar * (Qvalue*Pstar/P)^theta)^(1/theta) * P * Evalue;
%define two cobb douglas like coefficients
a = alphan * (theta - 1)/theta;
b = alphak * (theta - 1)/theta;
total = Evalue * M^((a + b)/(1-a-b)) * (a/w)^(a/(1-a-b)) * (b/r)^(b/(1-a-b));
domestic = (1/(X*xi^(1-theta)* (1 - tau)^theta * (Qvalue*Pstar/P)^theta *Cstar + 1))^((theta - 1)/theta) * P * total;
export = X  * Qvalue * Pstar * Cstar^(1/theta) * (1/(xi + xi^theta * (1 - tau)^(-theta) * ...
    (P/(Qvalue*Pstar))^theta/Cstar))^((theta - 1)/theta) * total;
end