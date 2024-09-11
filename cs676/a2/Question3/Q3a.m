% A2Q3(a)

% Initial input values 
sigma = 0.25;   % volatility
r = 0.03;       % risk-free interest rate
mu = 0.18;      % drift
T = 1;          % expiry time
S0 = 90;        % initial asset price
K = S0;         % strike price
N = 250;        % number of timesteps
dt = T/N;       % time interval
M = 10000;      % number of simulations

% at time t(N-1)
S = linspace(70,110,100);   % series of stock price at time t(N-1)
[C,P] = blsprice(S,K,r,dt,sigma);    % corresponding payoff at t(N-1)
V = C+P;

% simulate change of S
dS = zeros(M,length(S));
for i = 1:length(S)
    dS(:,i) = mu*S(i)*dt + sigma*S(i)*randn(M,1)*sqrt(dt);
end

% at time t(N)=T
ST = zeros(M,length(S));
for i = 1:length(S)
    ST(:,i) = S(i) + dS(:,i);
end
VT = max(ST-K,0) + max(K-ST,0);

% compute the delta hedging positions at time t(N-1) for each value of S
delta = zeros(length(S),1);
for i = 1:length(S)
    delta(i) = mean((VT(:,i) - V(i))./dS(:,i));
end

plot(S,delta)