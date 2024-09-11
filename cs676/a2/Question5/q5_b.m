

T = 1.00;               % time to expiry
sigma = 0.15;           % volatility
r = 0.03;               % risk free interest rate
S_init = 90;
K = S_init;
%jump size: log normal distribution
p_u = 0.4;             % probability of up jump
lambda = .1;             % jump size arrival rate lambda = 0.1



u_u = 0.32;             % parameter for up movement
u_d = 0.3;              % parameter for down movement
N_sim = 25000;
N=800;
delt = T/N;         % delt = 1/1000

jump_up_mean = u_u;
jump_down_mean = u_d;
%
% compenstaed drift E[J-1]
%

%kappa value that we found int question 5 (a)
kappa = p_u/(1-u_u) + (1-p_u)/(u_d + 1) - 1;
%
% compensated drift for X = log(S)

drift = r - 0.5*sigma^2 - lambda*kappa;

%
% X = log(S)
%

X_old(1:N_sim,1) = log(S_init);
X_new(1:N_sim,1) = zeros(N_sim,1);

jump_chek = zeros(N_sim,1);
jump_size = zeros(N_sim,1);
jump_mask = zeros(N_sim,1);

for i = 1:N %timestep loop

    
    jump_chek(:,1) = rand(N_sim,1);                         % first check to determine lambda*delt
    jump_chek2(:, 1) = rand(N_sim, 1);                      % second check to determine p_u
    jump_mask(:,1) = (jump_chek(:,1) <= lambda*delt);       
    jump_mask2(:, 1) = (jump_chek2(:, 1) <= p_u);           
    jump_mask3(:, 1) = (jump_chek2(:, 1) > p_u);            
    jump_size = (jump_mask2 .* exprnd(u_u, N_sim, 1)) -...  % determine the jump_size
        (jump_mask3 .* exprnd(u_d, N_sim, 1));
    jump_size = jump_size .* jump_mask;                     
    
    X_new(:,1) = X_old(:,1) + drift*delt + sigma*sqrt(delt)*randn(N_sim,1)+...
        jump_size(:,1);
    X_old(:,1) = X_new(:,1);
end %timestpe loop

S(:,1) = exp(X_new(:,1));

n_bin = 200;
hist(S,n_bin);

W(:,1) = max(K-S,0);
%European put option
Ep = exp(-r*T)*mean(W(:,1));
disp(sprintf('put option price is : %.5g\n', Ep));


