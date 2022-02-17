function [domestic,export] = sales(X, Qvalue, Evalue, theta, Cstar, alphan, alphak, w, r)
% Domestic and export sales for firm decision

% Why I create the function? 
% Creating multiple outputs in the profit function force you to create many
% extra variables to store the variables and main function will not be as aesthetic

M = (1 + X * (Qvalue^theta) * Cstar)^(1/theta) * Evalue;
n = (w * (w*alphak/(r * alphan))^(-alphak*(theta - 1)/theta) * ...
   theta/(alphan * M * (theta - 1))) ^ (theta/(alphan * (theta - 1) + alphak * (theta - 1) - theta));

k = w * alphak * n/(r * alphan);

% Calculate domestic and foreign sales
domestic = Evalue * ((1/(X * Qvalue^theta * Cstar + 1)) * n^alphan * k^alphak)^((theta - 1)/theta);
export =  X * Qvalue * Evalue * Cstar^(1/theta) * ((Qvalue^theta * Cstar/(Qvalue^theta * Cstar + 1)) * n^alphan * k^alphak)^((theta - 1)/theta);

end

