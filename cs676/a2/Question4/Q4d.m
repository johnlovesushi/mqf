

S_init = 99;    % initial stock price
K = 100;        % strike price
T = 1.00;       % time to expiry
sigma = 0.25;   % volatility
r = 0.03;       % risk free rate
mu = 0.18;      % expected stock return
N = [100,200,400,800];  % rebalancing times
delt = T ./ N;          % time interval
 
% calculating option price using MATLAB function blsprice
C0 = blsprice(S_init, K, r, T, sigma);

% B0 = C0 at t = 0
B0 = C0;
% simulating sample path
pie = zeros(80000,4);   % pie value initialization
for j = 1:4
    S_sim = zeros(80000,N(j)+1);
    B_sim = zeros(80000,N(j)+1);
    S_sim(:,1) = S_init;
    B_sim(:,1) = B0;
    for i = 1:N(j)
        % calculating  relative hedging error value
        % if currently holding the stock & next period S < K, then sell the stock
        % if currently not holding the stock & next period S > K, then buy stock
        S_sim(:,i+1) = S_sim(:,i).*exp((mu-0.5*sigma^2)*delt(j) + sigma*sqrt(delt(j))*randn(80000,1));
        B_sim(:,i+1) = B_sim(:,i).*exp(r*delt(j)) - ...
            (S_sim(:,i+1) > K & S_sim(:,i) <=K).*S_sim(:,i+1) +(S_sim(:,i+1) < K & S_sim(:,i) >= K).*S_sim(:,i+1);
    end
    pie(:,j) = exp(-r*T).*(-max(S_sim(:,N(j)+1)-K,0) + (S_sim(:,N(j)+1) > K) .*S_sim(:,N(j)+1) +B_sim(:,i))./C0;
end

%mean, standard deviation, VaR and CVaR initialization
mea = zeros(4,1);
sd = zeros(4,1);
VaR = zeros(4,1);
CVaR = zeros(4,1);

for j = 1:4
    mea(j) = mean(pie(:,j));
    sd(j) = std(pie(:,j));
    [VaR(j),CVaR(j)] = dVaRCVaR(pie(:,j),0.95); % call the function that we write in question 3

end
%generate a table 
Table = table(mea,sd,VaR,CVaR);
Table.Properties.VariableNames = {'mean','sd','VaR','CVaR'};
Table.Properties.RowNames = {'N=100','N=200','N=400','N=800'};

histogram(pie(:,4));
