function [ V0, S, Delta ] = binomialDeltaStraddle( S0, r, sigma, T, N, K )
% option pricing and delta computing
% V0 - intital option price
 
% S - a (N + 1) by (N + 1) matrix representing the underlying price on the
% binomial lattice nodes over equally spaced N time periods in [0,T].
 
% Delta - a N by N matrix representing the delta hedging positions on the
% binomial lattice nodes over equally spaced N time periods in [0,T].
 
% S0 - initial asset price
% K - strike
% T - expiry time
% r - interest rate, constant
% sigma - volatility, constant
% N - number of timesteps
 
 
% tree parameters?
 
delt = T / N;
u = exp(sigma * sqrt(delt) );
d = 1 ./u;
a = exp(r*delt);
p = (a - d)/(u - d);
 
 
S = zeros(N+1, N+1); % define a N+1 by N+1 matrix
 
% store prices into S
 
S(1, 1) = S0;
 
for i = 1 : N
    s = S0*u.^([i : -1 : 0]').*d.^([0 : i]');
    S(1:i+1, i+1) = s;
end
 
%
% payoff at t = T
%
W = max( s - K, 0) + max( K - s, 0);
 
% backward recursion to calculate delta on the binomial tree lattice
 
Delta = zeros(N, N);
Delta(1:N, N) = ( W(1:N) - W(2:N+1) ) ./ ((u-d) * S(1:N, N));
 
for i = N : -1 : 2
    W = exp(-r*delt)*( (1-p)*W(2:i+1) + p*W(1:i) );
    Delta(1:i-1, i-1) = ( W(1:i-1) - W(2:i) ) ./ ((u-d) * S(1:i-1, i-1));
end
 
% calculate option price
V0 = exp(-r*delt)*( (1-p)*W(2) + p*W(1) );
S;
Delta;
 
end
