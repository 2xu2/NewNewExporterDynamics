% Modified from the Julia code from the Quantecon package
% https://github.com/QuantEcon/QuantEcon.jl/blob/013ec2682794abccc5868392b7d7fe14b41f13dc/src/markov/markov_approx.jl#L21
function [s, Pi] = tauchen_quantecon(mu,rho,sig,N)
if N == 1
    s = 0;
    Pi = 1;
else
    % mu: Mean of AR(1) process
    % n_std: number of std each side should span
    n_std = 3;

    % Get discretized space
    a_bar = n_std * sqrt(sig^2/ (1 - rho^2));
    d = 2 * a_bar/(N - 1);
    y = -a_bar:d:a_bar;

    Pi = zeros(N, N);

    for row = 1:N
        % end points
        Pi(row, 1) = normcdf((y(1) - rho*y(row) + d/2) / sig);
        Pi(row, N) = 1 - normcdf((y(N) - rho*y(row) - d/2) / sig);

        % middle columns
        for col = 2:N-1
            Pi(row, col) = (normcdf((y(col) - rho*y(row) + d/2)/sig) - ...
                normcdf((y(col) - rho*y(row) - d/2)/sig));
        end
    end

    % center process around its mean (wbar / (1 - rho)) in new variable
    s = y + mu/(1 - rho);

    % renormalize
    Pi = Pi./sum(Pi,2);
end
end
