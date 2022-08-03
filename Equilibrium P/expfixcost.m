function f = expfixcost(X1, X2, model)
% Not currently used in the main code, could be used to extend the model
% X1 is state of firm whether to export at last period
% X2 is state of firm whether to export at current period
% model controls for identification for fixed cost
% model = 1: standard sunk cost
% model = 2; sunk cost with high elasticity
% model = 3; sunk cost with low elasticity
if model == 1
    f_E = 0.961362;
    f_C = 0.0472772;
elseif model == 2
    f_E = 0.736;
    f_C = 0.034;
elseif model == 3
    f_E = 1.312;
    f_C = 0.075;
end

f = (1-X1) * f_E * X2 + X1 * f_C * X2;
end

