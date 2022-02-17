function f = expfixcost(X1, X2, model)
% To be extended
% X1 is state of firm whether to export at last period
% X2 is state of firm whether to export at current period
% model controls for identification for fixed cost
% model = 1: standard sunk cost
% model = 2; sunk cost with high elasticity
% model = 3; sunk cost with low elasticity
if model == 1
    f0 = 0.961362;
    f1 = 0.0472772;
elseif model == 2
    f0 = 0.736;
    f1 = 0.034;
elseif model == 3
    f0 = 1.312;
    f1 = 0.075;
end

if X2 == 0
    f = 0;
elseif X1 == 0
    f = f0; % new exporter
else
    f = f1; % existing exporter
end
end

