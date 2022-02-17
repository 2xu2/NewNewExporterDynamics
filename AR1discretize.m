function [m,trans] = AR1discretize(N, mu, rho, esigma)
% Monte-Carlo Approximation for AR(1) process: y_t = (1 - \rho) * \mu + \rho * y_t-1 + \episilon_t
% Method from Tauchen (1986), modified from Kim Rhul (2013)
% N is the number of nodes for Z
% mu is unconditional mean for y
% rho correlation
% esigma std. of epsilon (mean of epsilon is 0)

% Calculate standard deviation of the AR(1) processes
x = mu * (1 - rho);
sigma = esigma/sqrt(1 - rho^2);

% Define the grids as in Adda and Cooper
% p = 1/(2*N):1/N:1-1/(2*N);
% m = norminv(p,mu,sigma);
spread = 8;
m = [];
m(1) = mu - spread * sigma;
for i = 2:N
    m(i) = m(i - 1) + 2 * spread * sigma/(N - 1);
end

% Compute transition probabilities (Prob(y_t = z_j|y_t-1 = z_i))
transition = [];
for i = 1:N
    transition(i,1) = normcdf(((m(1) - rho * m(i) - x + (m(2) - m(1))/2))/esigma);
    for j = 2:(N-1)
        transition(i, j) = normcdf((m(j) - rho * m(i) - x + (m(j + 1) - m(j))/2)/esigma) ...
        - normcdf((m(j) - rho * m(i) - x - (m(j) - m(j - 1))/2)/esigma);
    end
    transition(i, N) = 1 - normcdf((m(N) - rho * m(i) - x - (m(N) - m(N - 1)/2))/esigma);
end

% normalize row sum to 1
trans = [];
for i = 1:N
    trans(i,:) = transition(i,:)/sum(transition(i,:));
end
end

