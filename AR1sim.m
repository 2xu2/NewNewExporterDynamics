function y = AR1sim(T, Z, Zprob)
% Simulation for AR(1) process using tauchenHussey.m
% I follow the process of simulating AR(1) from Adda and Cooper
% T is the full time span, Z and Zprob from tauchenHussey.m
oldind = 1; % initial state (can be chosen arbitrarily of all states)
y = zeros(1, T);
for t = 1:T
    u = rand(1, 'double'); % u governs the random shock
    sum = 0; % sum is used to decide which state z is at
    ind = 1; % index of all period
    while sum < u
        sum = sum + Zprob(oldind, ind);
        ind = ind + 1;
    end
    y(t) = Z(ind - 1); % store the state in y
    oldind = ind - 1; % update initial state
end
end