

S_init = 90;            % initial stock price
K = 90;                 % strike price
N_sim = 10000;                      % number of simulation
T = 1.00;                           % time to expiry
sigma = 0.25;                       % volatility
r = 0.03;                           % risk free rate
mu = 0.18;                          % expected stock return
N = 250;                            % timesteps
delt = T / N;                       % time interval


%S_old = S_init * ones(N_sim, 1);    % a vector to store simulated stock price
%S_new = zeros(N_sim, 1);            % a vector to store simulated stock price
%S_sim = S_old;                      % a matrix to store all sample paths

% calculating option price using MATLAB function blsprice
[C0,P0] = blsprice(S_init, K, r, T, sigma);
V0 = C0 + P0;
% MC to generate the S prices by BS SDE formula
S_sim = zeros(N_sim,N+1);
S_sim(:,1) = S_init;
for i = 1:N
    S_sim(:,i+1) = S_sim(:,i).*exp((mu-0.5*sigma^2)*delt + sigma*sqrt(delt)*randn(N_sim,1));
end
VN = max(S_sim(:, N+1) - K, 0) + max(K - S_sim(:, N+1),0);    %calculate the payoff 

Delta_sim = zeros(N_sim,N);
for i = 1:N
    Delta_sim(:,i) = blsdelta(S_sim(:,i),K,r,delt*(N-i+1),sigma);
end

 %for no hedging

B0 = V0;                              %get initial bond value

pie_T_0 = -VN + B0 .* exp(r*T);         %get time T pie value
PL_0 = (pie_T_0) .* exp(-r*T) ./ V0; 
figure('Name','No Hedging');
hist(PL_0,50);


% hedging position rebalanced at n = 0
B0 = V0 - Delta_sim(:,1) .* S_sim(:,1); % computing B0
 % option value at expiry time
pie_T_1 = -VN + Delta_sim(:,1) .* S_sim(:,N+1) + B0 .* exp(r*T);
 % value of portfolio at expiry T
PL_1 = pie_T_1 .* exp(-r*T) ./ V0; % relative hedging error
% plot histogram
figure('Name',' Hedging at n=0');
hist(PL_1,50);




% daily hedging
B1 = B0; % initial hedging position

% computing first hedging
for n = 2 : 1 : N
    B1 = exp(r * delt) .* B1 + S_sim(:, n) .* (Delta_sim(:, n-1) - Delta_sim(:, n));
     % computing hegding positon of bonds
end
pie_T_2 = -VN + Delta_sim(:, N) .* S_sim(:, N+1) + B1*exp(delt*r);
PL_2 = pie_T_2 .* exp(-r*T) ./ V0;

figure('Name','Daily Hedging ');
hist(PL_2,50);



% hedging position rebalanced weekly
Bw = B0;

for n = 6 : 5 : N
    Bw = exp(r * delt * 5) .* Bw + S_sim(:, n) .* (Delta_sim(:, n-5)-Delta_sim(:, n));
     % computing hegding positon of bonds
end
% option value at expiry time
pie_T_3 = -VN + Delta_sim(:, n) .* S_sim(:, N+1) + Bw*exp(r*(N-n)*delt);
PL_3 = pie_T_3 .* exp(-r*T) ./ V0;
% plot histogram

figure('Name','Weekly Hedging ');
hist(PL_3,50);
 
 
% hedging position rebalanced monthly
Bm = B0;
for n = 21 : 20 : N
    Bm = exp(r * delt * 20) .* Bm + S_sim(:, n) .* (Delta_sim(:, n-20) - Delta_sim(:, n) );
     % computing hegding positon of bonds
end
% option value at expiry time
pie_T_4 = -VN + Delta_sim(:, n) .* S_sim(:, N+1) + Bm * exp(r * delt * (N-n));
PL_4 = pie_T_4 .* exp(-r*T) ./ V0;
% plot histogram

figure('Name','Monthly Hedging ');
hist(PL_4,50);







