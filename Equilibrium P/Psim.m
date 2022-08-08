function [P,Pstar, iterate] = Psim(xi_E, xi_C, rho_EE, rho_CC, f_E, f_C, Erho,Eesigma, EN, Qrho,Qesigma, QN, theta, Cstar, alphan, alphak, w, r, tau)
% Simulate for aggregate price level in domestic and foreign market
% Iterate to find equilibrium prices

%=========================================================================
% Initialize
%=========================================================================
T = 1000; % total time span to iterate over
plants = 2000; % number of plants to simulate 
Tdrop = 50; % Initial periods to drop
R = 1/(1+r);
% sales of a nonexporting median plant

% Initial guesses of the prices
MC = (w/alphan)^alphan * (r/alphak)^alphak;
P = theta/(theta - 1) * MC;
Pstar = theta/((theta - 1) * (1-tau)) * xi_E * MC;
% domestic_m will be used to calculate fixed cost
[domestic_m, ~] = sales(0, xi_E, 1, 1, theta, Cstar, alphan, alphak, w, r, 0, P, Pstar);

% Discretized AR(1) processes
[Estate_or, Eprob] = tauchen_quantecon(0, Erho,Eesigma, EN);
[Qstate_or, Qprob] = tauchen_quantecon(0, Qrho,Qesigma, QN);

% Note we have ln processes so we need to take exponentials
Estate = exp(Estate_or);
Qstate = exp(Qstate_or);

% joint transition probability for the random states
trans_joint = kron(Eprob, Qprob);

% Create grids
States = []; 
for j = 1:length(Estate)
    for k = 1:length(Qstate)
        % iterate over all possible states of two shocks
        States = [States; Estate(j), Qstate(k)];
    end
end

% Initialize value functions
% calculate profit in next period with different export states
V_NX = States; % Nonexporter
V_EX = States; % Exporter with high iceberg cost
V_CX = States; % Exporter with low iceberg cost
% Export states
V_NX(:,3) = 0; 
V_EX(:,3) = 1;
V_CX(:,3) = 1;
% Assign value of xi
V_NX(:,4) = 100; % Doesn't matter, pick any number over xi_C
V_EX(:,4) = xi_E;
V_CX(:,4) = xi_C;
% Calculate profits if not exporting next period
V_NX(:,5) = profit(0, 100, States(:,1), States(:,2), theta, Cstar, alphan, alphak, w, r, tau, P, Pstar);
V_EX(:,5) = profit(1, xi_E, States(:,1), States(:,2), theta, Cstar, alphan, alphak, w, r, tau, P, Pstar);
V_CX(:,5) = profit(1, xi_C, States(:,1), States(:,2), theta, Cstar, alphan, alphak, w, r, tau, P, Pstar);
% Calculate profits if exporting next period
V_NX(:,6) = profit(0, 100, States(:,1), States(:,2), theta, Cstar, alphan, alphak, w, r, tau, P, Pstar) - ...
    f_E * domestic_m;
V_EX(:,6) = profit(1, xi_E, States(:,1), States(:,2), theta, Cstar, alphan, alphak, w, r, tau, P, Pstar) - ...
    f_C * domestic_m;
V_CX(:,6) = profit(1, xi_C, States(:,1), States(:,2), theta, Cstar, alphan, alphak, w, r, tau, P, Pstar) - ...
    f_C * domestic_m;
% Guess of initial value functions
V_NX(:,7) = V_NX(:,5);
V_EX(:,7) = V_EX(:,5);
V_CX(:,7) = V_CX(:,5);
% Assign values for speed
V_NX(:,8:11) = 0;
V_EX(:,8:11) = 0;
V_CX(:,8:11) = 0;

%=========================================================================
% Now we loop until we find equilibrium prices
% Solve the model
%=========================================================================
iterate = 0;
error = 1;
P_last = 0;
Pstar_last = 0;
while error > 1e-5 && iterate < 100
    % Calculate profits if not exporting next period
    V_NX(:,5) = profit(0, 100, States(:,1), States(:,2), theta, Cstar, alphan, alphak, w, r, tau, P, Pstar);
    V_EX(:,5) = profit(1, xi_E, States(:,1), States(:,2), theta, Cstar, alphan, alphak, w, r, tau, P, Pstar);
    V_CX(:,5) = profit(1, xi_C, States(:,1), States(:,2), theta, Cstar, alphan, alphak, w, r, tau, P, Pstar);
    % Calculate profits if exporting next period
    V_NX(:,6) = profit(0, 100, States(:,1), States(:,2), theta, Cstar, alphan, alphak, w, r, tau, P, Pstar) - ...
        f_E * domestic_m;
    V_EX(:,6) = profit(1, xi_E, States(:,1), States(:,2), theta, Cstar, alphan, alphak, w, r, tau, P, Pstar) - ...
        f_C * domestic_m;
    V_CX(:,6) = profit(1, xi_C, States(:,1), States(:,2), theta, Cstar, alphan, alphak, w, r, tau, P, Pstar) - ...
        f_C * domestic_m;
    % Guess of initial value functions
    V_NX(:,7) = V_NX(:,5);
    V_EX(:,7) = V_EX(:,5);
    V_CX(:,7) = V_CX(:,5);

    % Value function iteration
    V = VFI_baseline(V_NX, V_EX, V_CX, rho_EE, rho_CC, trans_joint, R);

    %=========================================================================
    % Simulate price dynamics
    %=========================================================================
    % Shocks
    Esim = [];
    for i = 1:plants
        E = AR1sim(T, Estate, Eprob);
        Esim = [Esim; E];
    end
    Qsim = AR1sim(T, Qstate, Qprob);
    % initialize iceberg states for firms
    xi_sim = zeros(plants, T);
    xi_sim(:,1) = 100; % Assume everyone is a nonexporter at start, doesn't matter here
    % initialize export status for firms
    X_sim = zeros(plants, T);
    % initialize firm and aggregate price level
    domesticprice = zeros(plants, T);
    foreignprice = zeros(plants, T);
    aggdoprice = zeros(T,1);
    aggfoprice = zeros(T,1);

    % store export and iceberg cost status from last and current period
    last_X = zeros(plants, 1);
    last_xi = xi_sim(:,1);

    % first period
    for i = 1:plants
        [domesticprice(i,1), foreignprice(i,1)] = price(0, 100, Esim(i, 1), Qsim(1),...
            theta, alphan, alphak, w, r, tau);
    end
    foreignprice = foreignprice(:,1) * Qsim(1);
    
    % simulation
    for j = 2:T
        for i = 1:plants
            X_sim(i, j) = V((V(:, 1) == Esim(i, j-1)) & (V(:, 2) == Qsim(j-1)) & ...
                (V(:, 3) == last_X(i)) &  (V(:, 4) == last_xi(i)), 10);
            xi_sim(i, j) = icebergsim(X_sim(i, j), last_xi(i), xi_C, xi_E, rho_CC, rho_EE);
            [domesticprice(i,j), foreignprice(i,j)] = price(X_sim(i, j), xi_sim(i, j), Esim(i, j), Qsim(j), ...
                theta, alphan, alphak, w, r, tau);
        end
        aggdoprice(j) = (sum(domesticprice(:,j).^(1-theta)))^(1/(1-theta));
        aggfoprice(j) = (sum(nonzeros((foreignprice(:,j)) ./ Qsim(j)).^(1-theta)))^(1/(1-theta));
    end
    P = mean(aggdoprice((Tdrop+1):end));
    Pstar = mean(aggfoprice((Tdrop+1):end));
    error = max(abs(P - P_last), abs(Pstar - Pstar_last));
    P_last = P;
    Pstar_last = Pstar;
    % update domestic_m (?)
    [domestic_m, ~] = sales(0, xi_E, 1, 1, theta, Cstar, alphan, alphak, w, r, 0, P, Pstar);
    iterate = iterate + 1;
end
end