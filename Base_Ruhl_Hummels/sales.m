function [domestic,export] = sales(X, xi, Evalue, Qvalue, theta, Cstar, alphan, alphak, w, r)
% Domestic and export sales for firm decision

% Not fastest, but clear
M = (1 + X * xi^(1 - theta) * Cstar * (Qvalue^theta))^(1/theta) * Evalue;
%define two cobb douglas like coefficients
a = alphan * (theta - 1)/theta;
b = alphak * (theta - 1)/theta;
total = Evalue * M^((a + b)/(1-a-b)) * (a/w)^(a/(1-a-b)) * (b/r)^(b/(1-a-b));
domestic = (1/(X*xi^(1-theta)*Qvalue^theta*Cstar + 1))^((theta - 1)/theta) * total;
export = X  * Qvalue * Cstar^(1/theta) * (1/(xi + xi^theta * Qvalue^(-theta)/Cstar))^((theta - 1)/theta) * total;
end