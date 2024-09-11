randn('state',100);
%
T = 1.00;           %expiry time
sigma = 0.25;       %voliatility
mu = .018;          %P measure drift
S_init = 90;       %initial value
N_sim = 10000;      % number of simulations
N=250;              % number of timesteps
delt = T/N;         % timestep
r = 0.03;           % risk-free interest rate
K = S_init;

drift = r*delt;     %under risk neutrail
sigma_sqrt_delt = sigma*sqrt(delt);

u = exp(sigma*sqrt(delt));
d = 1/u;
S = linspace(70,110,100);    %the stock price at T = N-1


delta_S = 0.01;     %initialized the delt_S value
%delta_S(1,:) = S(1,:).*(drift + sigma_sqrt_delt*randn(1,100));

%get simulated V0 price
%V = zeros(1,100);           % straddle option payoff initialization

%V(1,:) = max(delta_S(1,:) + delta_S - K,0) +  max(K- delta_S - S(1,:),0);
%V0_tao = zeros(1,100);
%V0_tao(1,:) = exp(-r*delt*(N-1)).*V(1,:);

%get V0 price

V0 = zeros(1,100);
[call0,put0] = blsprice(S(1,:), K, r, delt, .25);
V0 = call0 + put0;

[call,put] = blsprice(S(1,:)+delta_S, K, r, delt, .25);
V = call +put;
delta = zeros(1,100);
delta = (V(1,:) - V0(1,:))./delta_S;
%ratio = zeros(1,100);
%ratio(1,:) = delta(1,:)./S(1,:);
plot(S,delta);
title('Figure. delta values versus S from $70 to $110')
xlabel('S');
ylabel('delta');

