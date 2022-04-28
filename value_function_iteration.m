function [outputArg1,outputArg2] = value_function_iteration(gamma, beta, alpha, delta, Xlow, Xhigh, step)
% This function performs valkue function iteration for dynamic programming

% List of X (total income)
Xlist = Xlow:step:Xhigh;
v0 = 0; % Initial value (v_0(x))
c_L = Xlist(1);
% Initialize list of value function
V = [];
for i = 1:length(Xlist)
    V(i, 1) = util_CRRA(Xlist(i), gamma);
end

% Value function iteration
for epoch = 1:1000
    for i = 1:length(Xlist)
        c_H = Xlist(i);
        Clist = c_L:step:c_H;
        for j = 1:length(Clist)
            nextX = delta*(c_H - Clist(j));
            nextV = 
        end
    end
end
        V(i, epoch + 1) = 

