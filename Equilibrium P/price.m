function [do_p, fo_p] = price(X, xi, Evalue,Qvalue, theta, alphan, alphak, w, r, tau)
% domestic and foreign price of a firm
% The formula requires alphan + alphak = 1
% foreign price is 0 if the firm is not exporting

MC = Evalue^(theta/(1-theta)) * (w/alphan)^alphan * (r/alphak)^alphak;
do_p = theta/(theta - 1) * MC;
fo_p = X * theta/((theta - 1) * (1-tau)) * xi * MC/Qvalue;

end 