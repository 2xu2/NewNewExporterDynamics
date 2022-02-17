%==========================================================================
% Simulation
% 1914 plants and 420 quarters in original paper
%==========================================================================
% Simulation for the AR(1) processes 
% Note we have ln processes so we need to take exponentials of the processes
% We have 1914 plants and 400 periods
% Note we have to annualize the results (taking average for each year)
% E vary for each plant and Q is identical for all plants
Esim = [];
for i = 1:plants
    E = AR1sim(T, Estate_or, Eprob);
    Esim = vertcat(Esim, exp(E)); % transform to ln process
end
Qsim = AR1sim(T, Qstate_or, Qprob);
Qsim = exp(Qsim); % transform to ln process

% initialize policy function for firms
policy = ones(plants, 1);
% initialize domestic and foreign sales
domestic = zeros(plants, T);
exports = zeros(plants, T);

for i = 1:plants
    for j = 2:T
        policy(i, j) = V((V(:, 1) == policy(i, j - 1)) & (V(:, 2) == Esim(i, j)) & ...
            (V(:, 3) == Qsim(j)), 9);
    end
end


% obtain domestic and foreign sales
for i = 1:plants
    for j = 1:T
        [domestic(i, j),export(i, j)] = sales(policy(i, j), Qsim(j), Esim(i, j),...
             theta, Cstar, alphan, alphak, w, r);
    end
end

% Remove first 20 periods
policy = policy(:,21:end);
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
meanexport = mean(exportsales(exportsales > 0));

% save firm's last period decision
last = policy(:,1);

% Now obtain starter and stopper ratio
current = [];
starterrate = [];
stopperrate = [];
for j = 2:(T - 20)
    current = policy(:,j);
    starter = 0;
    stopper = 0;
    totalnonexporter = 0;
    totalexporter = 0;
    for i = 1:plants
        if last(i) == 0
            totalnonexporter = totalnonexporter + 1;
            if current(i) == 1
                starter = starter + 1;
            end
        elseif last(i) == 1
            totalexporter = totalexporter + 1;
            if current(i) == 0
                stopper = stopper + 1;
            end
        end
    end
    starterrate(j) = starter/totalnonexporter;
    stopperrate(j) = stopper/totalexporter;
    last = current;      
end
            

% For starter and stopperrate equivalent to take total average
averagestarter = mean(starterrate);
averagestopper = mean(stopperrate);