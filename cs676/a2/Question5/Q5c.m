% A2Q5(c)

r = 0.03;       % risk-free interest rate
T = 1;          % expiry time
S0 = 90;        % initial asset price
K = linspace(70,120,20);    % strike prices

% cumput put option values from JumpPut function
P = JumpPut(800,250000,K);

% Find implied volatility from blsimpv
ImpVol = blsimpv(S0, K, r, T, P,0.5, 0, [], {'Put'});

% Plot volatility vs strike price
plot(K,ImpVol);
xlabel('Strike Price');
ylabel('Implied Volatility');