% A2Q5(b)

function P = JumpPut(N,M,K)
% N is the number of timesteps
% M is the number of simulations
% K is the strike price
% P is the output put option price

% Initial input values 
sigma = 0.15;   % volatility
r = 0.03;       % risk-free interest rate
T = 1;          % expiry time
S0 = 90;        % initial asset price
muu = 0.32;     % up drift
mud = 0.3;      % down drift
pu = 0.4;       % probability of up jump
pd = 1 - pu;    % probability of down jump
lambda = 0.1;   % intensity of Poisson Process
dt = T/N;       % time interval
kappa = pu/(1-muu) + pd/(1+mud) - 1;

% drift for Xt = log(St)
drift = (r - lambda * kappa - 0.5*sigma^2);

% Simulate Y = log(J) double exponential and compute jumpsize
X = log(S0) * ones(M,1);

% In each time step, simulate Y = log(J) which follows double exponential 
% distribution.
% To decide whether a jump happends, simulate an uniform(0,1) r.v. V, 
% if V <= lambda*dt, a jump occurs, otherwise no jump occurs.
% Then compute jumpsize as Y*I(V<=lambda*dt), and use Euler method to 
% simulate X =log(S).
for i = 1:N
    U = rand(M,1);
    ind = (U <= pu);
    Y = ind.*exprnd(1/muu,M,1) - (1-ind).*exprnd(1/mud,M,1);
    V = rand(M,1);
    jump_ind = (V <= lambda*dt);
    jumpsize = jump_ind .* Y;
    X = X + drift*dt +  sigma*sqrt(dt)*randn(M,1) + jumpsize;
end

% Recover S = exp(X)
S = exp(X);

% Compute payoff at expiry T
payoff = max(K-S,0);

% Estimate option value
P = mean(exp(-r*T)*payoff);

end
