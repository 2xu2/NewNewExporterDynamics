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
    EN = 100;
	% QN is the grids (states) of markov chain to approximate AR(1) process of Q
    QN = 1;
    % Coefficients for iceberg cost
    % xi_E and xi_C are iceberg costs for exporting firms
    % rho_EE, rho_CC determines transition matrix of the markov chain
    xi_E = 1.6;
    xi_C = 1.2;
    rho_EE = 0.92;
    rho_CC = 0.92;
    % export fixed cost, f_E for new entry, f_C for continuing to export
    f_E = 0.9;
    f_C = 0.2;
    % theta is the elasticity of substitution
    theta = 5;
    % Cstar is the size of world aggregate demand relative to domestic demand
    Cstar = 0.146496;
    % alphan and alphak are parameters of the Copp-Douglas production function
    alphan = 0.45;
    alphak = 0.55;
    % tau is ad-valorem tariff
    tau = 0;
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
T = 100; % total time span to iterate over
plants = 2000; % number of plants to simulate
R = 1/(1 + r); % discount factor for value function

% sales of a nonexporting median plant
% domestic_m will be used to calculate fixed cost
[domestic_m,export_m] = sales(0, xi_E, 1, 1, theta, Cstar, alphan, alphak, w, r, 0, 1, 1);

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

%=========================================================================
% P and P_star
% find equilibrium aggregate price levels for both markets
%==========================================================================
[P,Pstar, iterate, Ppath] = Psim(xi_E, xi_C, rho_EE, rho_CC, f_E, f_C, Erho,Eesigma, ...
    EN, Qrho,Qesigma, QN, theta, Cstar, alphan, alphak, w, r, tau);
toc

tic
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
% We have 1914 plants and 400 periods
% Note we have to annualize the results (taking average for each year)
% E vary for each plant and Q is identical for all plants
tic
rng(2); % Set random seed
result = zeros(6, 1); 
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
% initialize prices
domesticprice = zeros(plants, T);
foreignprice = zeros(plants, T);
priceratio = zeros(T, 1);
avgforeignprice = zeros(T, 1);
% store export and iceberg cost status from last and current period
last_X = zeros(plants, 1);
last_xi = xi_sim(:,1);
totalexporterlast = 0;
current = zeros(plants, 1);
starterrate = zeros(T, 1);
stopperrate = zeros(T, 1);
churning = zeros(T, 1);
exportershare = zeros(T, 1);
% first period
for i = 1:plants
    [domestic(i, 1),export(i, 1)] = sales(0, 100, Esim(i, 1), Qsim(1), ...
             theta, Cstar, alphan, alphak, w, r, tau, P, Pstar);
    [domesticprice(i,1), foreignprice(i,1)] = price(0, 100, Esim(i, 1),  ...
            theta, alphan, alphak, w, r, tau);
end

% simulation
for j = 2:T
    starter = 0;
    stopper = 0;
    totalexporter = 0;
    for i = 1:plants
        X_sim(i, j) = V((V(:, 1) == Esim(i, j-1)) & (V(:, 2) == Qsim(j-1)) & ...
             (V(:, 3) == last_X(i)) &  (V(:, 4) == last_xi(i)), 10);
        xi_sim(i, j) = icebergsim(X_sim(i, j), last_xi(i), xi_C, xi_E, rho_CC, rho_EE);
        [domestic(i, j),export(i, j)] = sales(X_sim(i, j), xi_sim(i, j), Esim(i, j), Qsim(j), ...
             theta, Cstar, alphan, alphak, w, r, tau, P, Pstar);
        exportsales(i, j) = export(i, j)/(domestic(i, j) + export(i, j));
        [domesticprice(i,j), foreignprice(i,j)] = price(X_sim(i, j), xi_sim(i, j), Esim(i, j), ...
            theta, alphan, alphak, w, r, tau);
        if (last_X(i) == 0) && (X_sim(i, j) == 1)
                starter = starter + 1;
        end
        if (last_X(i) == 1) && (X_sim(i, j) == 0)
                stopper = stopper + 1;
        end
        totalexporter = totalexporter + X_sim(i, j);
    end
    avgforeignprice(j) = mean(nonzeros(foreignprice(:,j) .* Qsim(j)));
    exportershare(j) = totalexporter/plants;
    starterrate(j) = starter/(plants - totalexporterlast);
    stopperrate(j) = stopper/totalexporterlast;
    churning(j) = 2 * (starter + stopper)/(totalexporterlast + totalexporter);
    last_X = X_sim(:, j);
    last_xi = xi_sim(:, j);
    totalexporterlast = totalexporter;
end
avgdomesticprice = mean(domesticprice, 1, 'omitnan')';


% generated moments
result(1) = mean(starterrate(2:end), 'omitnan');
result(2) = mean(stopperrate(2:end), 'omitnan');
result(3) = mean(exportsales(exportsales > 0));
result(4) = mean(churning(2:end), 'omitnan');
result(5) = mean(exportershare(2:end));
result(6) = mean(avgforeignprice./avgdomesticprice, 'omitnan');

toc
disp(result);

% Create plots
t = tiledlayout(6,1);
ax1 = nexttile;
plot(ax1, starterrate)
title(ax1,'starter rate')
ax2 = nexttile;
plot(ax2, stopperrate)
title(ax2,'stopper rate')
ax3 = nexttile;
plot(ax3, churning)
title(ax3,'churning rate')
ax4 = nexttile;
plot(ax4, exportershare)
title(ax4,'exporter shares')
ax5 = nexttile;
plot(ax5, avgdomesticprice)
title(ax5,'average domestic price')
ax6 = nexttile;
plot(ax6, avgforeignprice)
title(ax6,'average foreign price (in domestic currency)')
% Add shared title and axis labels
xlabel(t,'time periods')
