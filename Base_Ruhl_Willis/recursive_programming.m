clear;
rng(2); % Set random seed
%==========================================================================
% Initialization
%==========================================================================
parms=input('Specify set of parameters: self input=0, baseline sunkcost=1, high elasticity=2 ');
if parms==1
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
elseif parms == 2
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
T = 420; % total time span to iterate over (first four seasons are dropped)
plants = 1914; % number of plants to simulate
R = (1/(1 + r))^0.25; % discount factor for value function
Xspace = [0, 1]; % Define state spaces for export choice

% sales of a nonexporting median plant
% domestic_m will be used to calculate fixed cost
[domestic_m,export_m] = sales(0, 1, 1, theta, Cstar, alphan, alphak, w, r);

% Discretized AR(1) processes
% Different Methodologies listed below:
% [Estate, Eprob] = AR1discretize(EN,(-Eesigma^2/(2 * (1 - Erho^2))),Erho,Eesigma);
% [Qstate, Qprob] = AR1discretize(QN,(-Qesigma^2/(2 * (1 - Qrho^2))),Qrho,Qesigma);

% [Estate_or, Eprob] = AR1discretize(EN,0,Erho,Eesigma);
% [Qstate_or, Qprob] = AR1discretize(QN,0,Qrho,Qesigma);

%[Estate_or, Eprob] = mytauchen(0, Erho,Eesigma, EN);
%[Qstate_or, Qprob] = mytauchen(0, Qrho,Qesigma, QN);

[Estate_or, Eprob] = tauchen_quantecon(0, Erho,Eesigma, EN);
[Qstate_or, Qprob] = tauchen_quantecon(0, Qrho,Qesigma, QN);

% Note we have ln processes so we need to take exponentials of the processes
Estate = exp(Estate_or);
Qstate = exp(Qstate_or);

% joint transition probability for the random states
trans_joint = kron(Eprob, Qprob);

% Initialize value for value functions
State = []; % Stores states and profits of the states
for i = 1:length(Xspace)
    for j = 1:length(Estate)
        for k = 1:length(Qstate)
            % iterate over all possible states
            % calculate profit in next period of exporting and not exporting
            % make a guess on the value function with current profit (V(:,6))
            State = [State; Xspace(i), Estate(j), Qstate(k), ...
                profit(0, Estate(j), Qstate(k), theta, Cstar, alphan, alphak, w, r), ...
                profit(1, Estate(j), Qstate(k), theta, Cstar, alphan, alphak, w, r) - expfixcost(Xspace(i), 1, parms) * domestic_m];
        end
    end
end
loops = length(Xspace) * length(Estate) * length(Qstate); % size of states

%=========================================================================
% Value function iteration
% Here we need two values for choosing to or to not export
%==========================================================================
iterate = 0;
error = 1;
V_EX = zeros(loops, 1); % Stores the value function for exporters
V_NX = zeros(loops, 1); % Stores the value function for non-exporters
Policy = zeros(loops, 1); % Store the policy function (0, 1)
V_old = zeros(loops, 1); 
V_old(:,1) = State(:, 4); % Take initial guess of the value function
V_new = zeros(loops, 1); 
while error > 1e-5
    V_EX(1:loops/2) = State(1:(loops/2), 4) + R * (trans_joint * V_old(1:(loops/2)));
    V_NX(1:loops/2) = State(1:(loops/2), 5) + R * (trans_joint * V_old((loops/2+1):end));
    V_EX((loops/2 + 1):end) = State((loops/2 + 1):end, 4) + R * (trans_joint * V_old(1:(loops/2)));
    V_NX((loops/2 + 1):end) = State((loops/2 + 1):end, 5) + R * (trans_joint * V_old((loops/2+1):end));
    % Alternatively you can create a sparse diagonal matrix of trans_joint
    % by blkdiag(trans_joint, trans_joint) and multiply in one go, the 
    % method is slightly slower (2-second difference)
    % So here I just assign the two halves seperately. 
    % Determine policy function and maximize value function
    for j = 1:loops
        if V_EX(j) > V_NX(j)
            Policy(j) = 0;
            V_new(j) = V_EX(j);
        else
            Policy(j) = 1;
            V_new(j) = V_NX(j);
        end
    end
    error = max(abs(V_old - V_new));
    V_old(:, 1) = V_new(:, 1); % update value function
    iterate = iterate + 1;
    if iterate > 2000
        fprintf("Fail to converge in 2000 iterations")
        break
    end
end
fprintf("Converged in %d iterations.\n", iterate);
toc



