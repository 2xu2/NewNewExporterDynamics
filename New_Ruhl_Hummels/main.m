clear;
rng(2); % Set random seed
%==========================================================================
% Initialization
%==========================================================================
model=input('Specify model: self input=0, baseline sunkcost=1, high elasticity=2');
if model==1
    % r is the average observed interest rate
    r = 0.109;
    % Qrho and Qesigma defines the AR(1) process for lnQ
    Qrho = 0.826; 
    Qesigma = 0.036;
    % Erho and Eesigma defines the AR(1) process for lnEpsilon
    Erho = 0.872524;
    Eesigma = 0.115886;
    % EN is the grids (states) of markov chain to approximate AR(1) process of E
    EN = 455;
	% QN is the grids (states) of markov chain to approximate AR(1) process of Q
    QN = 11;
    % Coefficients for iceberg cost
    % xi_E and xi_C are iceberg costs for exporting firms
    % rho_EE, rho_CC determines transition matrix of the markov chain
    xi_E = 1.2;
    xi_C = 1.07;
    rho_EE = 0.92;
    rho_CC = 0.92;
    % export fixed cost, f_E for new entry, f_C for continuing to export
    f_E = 0.1;
    f_C = 0.01;
    % theta is the elasticity of substitution
    theta = 5;
    % Cstar is the size of world aggregate demand relative to domestic demand
    Cstar = 0.146496;
    % alphan and alphak are parameters of the Copp-Douglas production function
    alphan = 0.45;
    alphak = 0.55;
    %Set w so that median firm has 63 employees (nstable = 63)
    nstable = 63; % stable level of employment
    w = alphan *((theta-1)/theta) * nstable^(alphan *(theta-1)/theta-1); 
elseif model == 2
    % r is the average observed interest rate
    r = 0.109;
    % Qrho and Qesigma defines the AR(1) process for lnQ
    Qrho = 0.826; 
    Qesigma = 0.036;
    % Erho and Eesigma defines the AR(1) process for lnEpsilon
    Erho = 0.873;
    Eesigma = 0.087;
    % EN is the grids (states) of markov chain to approximate AR(1) process of E
    EN = 455;
	% QN is the grids (states) of markov chain to approximate AR(1) process of Q
    QN = 11;
    % theta is the elasticity of substitution
    theta = 7;
    % Cstar is the size of world aggregate demand relative to domestic demand
    Cstar = 0.135;
    % alphan and alphak are parameters of the Copp-Douglas production function
    alphan = 0.45;
    alphak = 0.55;
    %Set w so that median firm has 63 employees (nstable = 63)
    nstable = 63; % stable level of employment
    w = alphan *((theta-1)/theta) * nstable^(alphan *(theta-1)/theta-1);
else
    r = input('r: ');
    Qrho = input('Qrho: ');
    Qesigma = input('Qesigma: ');
    Erho = input('Erho: ');
    Eesigma = input('Eesigma: ');
    EN = input('EN: ');
    QN = input('QN: ');
    QN = input('theta: ');
    Cstar = input('Cstar: ');
    alphan = input('alphan: ');
    alphak = input('alphak: ');
    nstable = 63;
    w = alphan *((theta-1)/theta) * nstable^(alphan *(theta-1)/theta-1);
end

tic
%==========================================================================
% Preparations
%==========================================================================
% Some other parameters that don't vary with models
T = 20; % total time span to iterate over (first four seasons are dropped)
plants = 150; % number of plants to simulate
R = 1/(1 + r); % discount factor for value function
Xspace = [0, 1]; % Define state spaces for export choice

% sales of a nonexporting median plant
% domestic_m will be used to calculate fixed cost
[domestic_m,export_m] = sales(0, xi_E, 1, 1, theta, Cstar, alphan, alphak, w, r);

% Discretized AR(1) processes
[Estate_or, Eprob] = tauchen_quantecon(0, Erho,Eesigma, EN);
[Qstate_or, Qprob] = tauchen_quantecon(0, Qrho,Qesigma, QN);

% Note we have ln processes so we need to take exponentials of the processes
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
V_NX(:,5) = profit(0, 100, States(:,1), States(:,2), theta, Cstar, alphan, alphak, w, r);
V_EX(:,5) = profit(1, xi_E, States(:,1), States(:,2), theta, Cstar, alphan, alphak, w, r);
V_CX(:,5) = profit(1, xi_C, States(:,1), States(:,2), theta, Cstar, alphan, alphak, w, r);
% Calculate profits if exporting next period
V_NX(:,6) = profit(0, 100, States(:,1), States(:,2), theta, Cstar, alphan, alphak, w, r) - ...
    f_E * domestic_m;
V_EX(:,6) = profit(1, xi_E, States(:,1), States(:,2), theta, Cstar, alphan, alphak, w, r) - ...
    f_C * domestic_m;
V_CX(:,6) = profit(1, xi_C, States(:,1), States(:,2), theta, Cstar, alphan, alphak, w, r) - ...
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
% Value function iteration
% V(j, 1-4) are states we defined
% V(j, 5) is the profit this period with choosing to export next period
% V(j, 6) is the profit with choosing not to export next period
% V(j, 7) is the value function from last iteration
% V(j, 8) is the new value function choosing not to export next period
% V(j, 9) is the new value function choosing to export next period
% V(j, 10) is the policy function
% V(j, 11) is the maximized value function
%==========================================================================
V = VFI_baseline(V_NX, V_EX, V_CX, rho_EE, rho_CC, trans_joint, R);
toc

%==========================================================================
% Simulation
%==========================================================================
% Simulation for the AR(1) processes 
% Note we have ln processes so we need to take exponentials of the processes
% We have 1914 plants and 400 periods
% Note we have to annualize the results (taking average for each year)
% E vary for each plant and Q is identical for all plants
tic
rng(2); % Set random seed
result = zeros(3, 1); 
Esim = [];
for i = 1:plants
    E = AR1sim(T, Estate, Eprob);
    Esim = [Esim; E];
end
Qsim = AR1sim(T, Qstate, Qprob);
% initialize iceberg states for firms
xi_sim = zeros(plants, T);
xi_sim(:,1) = 100; % everyone is a nonexporter at start
% initialize export status for firms
X_sim = zeros(plants, T);
% initialize domestic and foreign sales
domestic = zeros(plants, T);
export = zeros(plants, T);
% initialize export-sales ratio of exporting plants
exportsales = zeros(plants, T);
% store export and iceberg cost status from last and current period
last_X = zeros(plants, 1);
last_xi = xi_sim(:,1);
current = zeros(plants, 1);
starterrate = zeros(T, 1);
stopperrate = zeros(T, 1);
% First period
for i = 1:plants
    [domestic(i, 1),export(i, 1)] = sales(0, 100, Esim(i, 1), Qsim(1), ...
             theta, Cstar, alphan, alphak, w, r);
    exportsales(i, 1) = export(i, 1)/(domestic(i, 1) + export(i, 1));
end

% simulate export choice and domestic and foreign sales
for j = 2:T
    starter = 0;
    stopper = 0;
    totalnonexporter = 0;
    for i = 1:plants
        X_sim(i, j) = V((V(:, 1) == Esim(i, j-1)) & (V(:, 2) == Qsim(j-1)) & ...
             (V(:, 3) == last_X(i)) &  (V(:, 4) == last_xi(i)), 10);
        xi_sim(i, j) = icebergsim(X_sim(i, j), last_xi(i), xi_C, xi_E, rho_CC, rho_EE);
        [domestic(i, j),export(i, j)] = sales(X_sim(i, j), xi_sim(i, j), Esim(i, j), Qsim(j), ...
             theta, Cstar, alphan, alphak, w, r);
        exportsales(i, j) = export(i, j)/(domestic(i, j) + export(i, j));
        % store current export status
        current(i) = X_sim(i,j);
        if last_X(i) == 0
            totalnonexporter = totalnonexporter + 1;
            if current(i) == 1
                starter = starter + 1;
            end
        else
            if current(i) == 0
                stopper = stopper + 1;
            end
        end
    end
    starterrate(j) = starter/totalnonexporter;
    stopperrate(j) = stopper/(plants - totalnonexporter);
    last_X = X_sim(:, j);
    last_xi = xi_sim(:, j);
end

% mean export-sales ratio
result(1) = mean(starterrate, 'omitnan');
result(2) = mean(stopperrate, 'omitnan');
result(3) = mean(exportsales(exportsales > 0));

toc
disp(result);



