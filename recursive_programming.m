clear;
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
    %w from paper
        w = (((theta-1)/theta)*alphan * (r * alphan/alphak)^(alphak *...
            (theta-1)/theta) * nstable^(((alphan + alphak)*(theta - 1) - ...
            theta)/theta))^(theta/(theta + alphak* theta - alphak)); % w I derived
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
    %w = alphan *((theta-1)/theta) * nstable^(alphan *(theta-1)/theta-1);
    w = (((theta-1)/theta)*alphan * (r * alphan/alphak)^(alphak *...
        (theta-1)/theta) * nstable^(((alphan + alphak)*(theta - 1) - ...
        theta)/theta))^(theta/(theta + alphak* theta - alphak));    
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
% [Estate, Eprob] = AR1discretize(EN,(-Eesigma^2/(2 * (1 - Erho^2))),Erho,Eesigma);
% [Qstate, Qprob] = AR1discretize(QN,(-Qesigma^2/(2 * (1 - Qrho^2))),Qrho,Qesigma);

% [Estate_or, Eprob] = AR1discretize(EN,0,Erho,Eesigma);
% [Qstate_or, Qprob] = AR1discretize(QN,0,Qrho,Qesigma);

[Estate_or, Eprob] = mytauchen(0, Erho,Eesigma, EN);
[Qstate_or, Qprob] = mytauchen(0, Qrho,Qesigma, QN);

% Note we have ln processes so we need to take exponentials of the processes
Estate = exp(Estate_or);
Qstate = exp(Qstate_or);
% % create grids of the two random states
% shockstate = [];
% for j = 1:length(Estate)
%     for k = 1:length(Qstate)
%         shockstate((j-1) * length(Qstate) + k, 1) = Estate(j);
%         shockstate((j-1) * length(Qstate) + k, 2) = Qstate(k);
%     end
% end
% joint transition probability for the random states
trans_joint = kron(Eprob, Qprob);

% Initialize value for value functions
V = [];
for i = 1:length(Xspace)
    for j = 1:length(Estate)
        for k = 1:length(Qstate)
            % iterate over all possible states
            % calculate profit in next period of exporting and not exporting
            % make a guess on the value function with current profit
            V = vertcat(V, [Xspace(i), Estate(j), Qstate(k), ...
                profit(0, Qstate(k), Estate(j), theta, Cstar, alphan, alphak, w, r), ...
                profit(1, Qstate(k), Estate(j), theta, Cstar, alphan, alphak, w, r) - expfixcost(Xspace(i), 1, parms) * domestic_m, ...
                profit(0, Qstate(k), Estate(j), theta, Cstar, alphan, alphak, w, r), ...
                0, 0, 0, 0]); % Pre-allocate some space
        end
    end
end

%=========================================================================
% Value function iteration
% Here we need two values for choosing to or to not export
%==========================================================================
for t = 1:100
    for j = 1:length(V)
        % V(j, 4) is the profit without export calculated
        % V(j, 5) is the profit with export calculated
        % V(j, 6) is the value function from last iteration
        % V(j, 7) is the new value function without export
        % V(j, 8) is the new value function with export
        % V(j, 9) is the policy function
        % V(j, 10) is the maximized value function
        if j <= length(V)/2
            V(j, 7) = V(j, 4) + R * (trans_joint(j,:) * V(1:(length(V)/2), 6));
            V(j, 8) = V(j, 5) + R * (trans_joint(j,:) * V((length(V)/2 + 1):length(V), 6));
        else 
            V(j, 7) = V(j, 4) + R * (trans_joint(j - length(V)/2,:) * V(1:(length(V)/2), 6));
            V(j, 8) = V(j, 5) + R * (trans_joint(j - length(V)/2,:) * V((length(V)/2 + 1):length(V), 6));
        end
        % Determine policy function and maximize value function
        if V(j, 7) > V(j, 8)
            V(j, 9) = 0;
            V(j, 10) = V(j, 7);
        else
            V(j, 9) = 1;
            V(j, 10) = V(j, 8);
        end
        V(j, 6) = V(j, 10); % update value function
    end
end



