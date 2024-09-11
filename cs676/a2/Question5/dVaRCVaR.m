% A2Q3(c)

function [VaR,CVaR] = dVaRCVaR(PL,beta)

M = length(PL);  

% sort P&L in increasing order
PL = sort(PL);

% find the index of M simulated P&L for beta (s.t. index/M <= 1-beta)
ind = floor(M*(1 - beta));

% compute VaR
VaR = PL(ind);

% compute CVaR by taking the mean of (P&L <= VaR)
CVaR = mean(PL(1:ind));

end