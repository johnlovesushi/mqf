randn('state',100);
%
T = 1.00;           %expiry time
sigma = 0.25;       %voliatility
mu = .018;          %P measure drift
S_init = 90;       %initial value
%N_sim = 10000;      % number of simulations
N=250;              % number of timesteps
delt = T/N;         % timestep
r = 0.03;           % risk-free interest rate
K = S_init;
drift = mu*delt;
sigma_sqrt_delt = sigma*sqrt(delt);


u = sigma_sqrt_delt;
d = 1/u;

q = (exp(r*delt) - d)/(u-d);
S_old = zeros(N,1);
S_new = zeros(N,1);
V = zeros(N+1,N+1);  %European straddle
S = zeros(N+1,N+1);

S(1:N,1) = S_init;
for i=1:N     %timestep loop
    % now, for each timestep, generate info for
    % all simulations
    
    S(:,i+1) = S(:,i) + ...
        S(:,i).*(drift + sigma_sqrt_delt*randn(N+1,1));
    
    S(:,i+1) = max(0.0,S(:,i+1));
        % check to make sure that S_new cannot be <0
    %S_old(:,1) = S_new(:,1);
        %
        %end of generation of all data for all simulations
        % for this timestep
end %timestep loop



V(:,N+1) = max(S(:,N+1) -K,0) + max(K-S(:,N+1),0);
%backward recursion
for i = N+1:-1:2
    V(1:i-1,i-1) = exp(-r*delt)*(q*V(1:i-1,i)+ (1-q)*V(2:i,i));
end


%measure delta

delta = zeros(N,N);
for i = N+1:-1:2
    delta(1:i-1,i-1) = (V(1:i-1,i) - V(2:i,i))./((u-d)*S(1:i-1,i-1));
end
%n_bin = 200;
%histogram(S_new,n_bin);
%histogram(V,n_bin);
%stndrd_dev = std(S_new);
%disp(sprintf('standard deviation: %.5g\n', stndrd_dev));
Sq4 = linspace(70, 110, 100);
plot(Sq4, delta); 
%mean_S = mean(S_new);
%disp(sprintf('mean: %.5g',mean_S));