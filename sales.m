function [domestic,export] = sales(X, Evalue, Qvalue, theta, Cstar, alphan, alphak, w, r)
% Domestic and export sales for firm decision

% Not fastest, but clear
M = (1 + X * (Qvalue^theta) * Cstar)^(1/theta) * Evalue;
%define two cobb douglas like coefficients
a = alphan * (theta - 1)/theta;
b = alphak * (theta - 1)/theta;
Y = Evalue * M^((a + b)/(1-a-b)) * (a/w)^(a/(1-a-b)) * (b/r)^(b/(1-a-b));
domestic = (1/(X*Qvalue^theta*Cstar + 1))^((theta - 1)/theta) * Y;
export = X  * Qvalue * Cstar^(1/theta) * (1/(1 + Qvalue^(-theta)/Cstar))^((theta - 1)/theta) * Y;
end