%==========================================================================
% Simulation
% 1914 plants and 420 quarters in original paper
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

% initialize policy function for firms
policy_sim = zeros(plants, T);
% initialize domestic and foreign sales
domestic = zeros(plants, T);
export = zeros(plants, T);

% simulate export choice and domestic and foreign sales
for i = 1:plants
    for j = 2:T
        policy_sim(i, j) = Policy((State(:, 1) == policy_sim(i, j - 1)) & (State(:, 2) == Esim(i, j)) & ...
            (State(:, 3) == Qsim(j)));
        [domestic(i, j),export(i, j)] = sales(policy_sim(i, j), Esim(i, j), Qsim(j), ...
             theta, Cstar, alphan, alphak, w, r);
    end
end


last = policy_sim(:,20); %Save the policy of the last removed period

% Remove first 20 periods
policy_sim = policy_sim(:,21:end);
domestic = domestic(:,21:end);
export = export(:,21:end);

% Annualize sales
domestic_y = zeros(plants, (T- 20)/4);
export_y = zeros(plants, (T- 20)/4);
for i = 1:(T - 20)/4
    for j = 1:4
        domestic_y(:, i) = domestic_y(:, i) + domestic(:,4*i -4 + j);
        export_y(:, i) = export_y(:, i) + export(:, 4*i -4 + j);
    end
end
% calculate the mean of export sales for each period
% initialize export-sales ratio of exporting plants
exportsales = zeros(plants, (T - 20)/4);
for i = 1:plants
    for j = 1:(T - 20)/4
    exportsales(i, j) = export_y(i, j)/(domestic_y(i, j) + export_y(i, j));
    end
end

% mean variance of domestic sales (plant size)
vardomestic = mean(var(domestic));

% mean export-sales ratio
result(3) = mean(exportsales(exportsales > 0));

% Now obtain starter and stopper ratio
starterrate = zeros(T-20, 1);
stopperrate = zeros(T-20, 1);
for j = 1:(T - 20)
    current = policy_sim(:,j);
    starter = 0;
    stopper = 0;
    totalnonexporter = 0;
    for i = 1:plants
        if last(i) == 0
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
    last = current;      
end
            
result(1) = mean(starterrate);
result(2) = mean(stopperrate);
toc
disp(result);
