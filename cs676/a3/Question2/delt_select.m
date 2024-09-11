function [result,iteration] = delt_select(S,N)
%sigma = 0.15;
alpha_parameter =0.8;    % constant related to alpha
r = .02;                 % risk free rate
T = 1;                   % time to expiry
K = 10;                  % strike price
S0 = 10;                 % initial stock price
%N = 25;
delt = T/N;              % time interval
%S = S4;
Large = 1e6;             % penalty coefficient
tolerance = 1/Large;     % toleranceeranceerance term
dnorm = 0.1;

%S = [0:0.1*K:0.4*K,...  %input S value
%    0.45*K:0.05*K:0.8*K,...
%    0.82*K:0.02*K:0.9*K,...
%    0.91*K:0.01*K:1.1*K,...
%    1.12*K:0.02*K:1.2*K,...
%    1.25*K:.05*K:1.6*K,...

%    1.7*K:0.1*K:2*K,...
%    2.2*K, 2.4*K, 2.8*K,...
%    3.6*K, 5*K, 7.5*K, 10*K];
%%%

Vs = max(K-S.^2, S.^2 - K)';
%Vs = max(K-S, S - K)';
m = length(S);  %number of grids in row

alpha_central = zeros(m,1);
beta_central = zeros(m,1);
alpha_forward = zeros(m,1);
beta_forward  = zeros(m,1);
alpha_backward = zeros(m,1);
beta_backward = zeros(m,1);
alpha = zeros(m,1);
beta = zeros(m,1);


for i = 2: m - 1
    %%central alpha and beta formula
    alpha_central(i) = alpha_parameter^2*S(i)/((S(i) - S(i-1))*(S(i+1) - S(i-1)))...
        -r*S(i)/(S(i+1) - S(i-1));
    beta_central(i) = alpha_parameter^2*S(i)/((S(i+1) - S(i))*(S(i+1) - S(i-1)))...
        +r*S(i)/(S(i+1) - S(i-1));
    %alpha_central(i) = sigma^2*S(i).^2/((S(i) - S(i-1))*(S(i+1) - S(i-1)))...
    %    -r*S(i)/(S(i+1) - S(i-1));
    %beta_central(i) = sigma^2*S(i).^2/((S(i+1) - S(i))*(S(i+1) - S(i-1)))...
    %    +r*S(i)/(S(i+1) - S(i-1));
    %%forward alpha and beta formula
    alpha_forward(i) = alpha_parameter^2*S(i)/((S(i) - S(i-1))*(S(i+1) - S(i-1)));
    beta_forward(i) = alpha_parameter^2*S(i)/((S(i+1) - S(i))*(S(i+1) - S(i-1)))...
        + r*S(i)/(S(i+1) - S(i));  % same as beta_central
    %alpha_forward(i) = sigma^2*S(i).^2/((S(i) - S(i-1))*(S(i+1) - S(i-1)));
    %beta_forward(i) = sigma^2*S(i).^2/((S(i+1) - S(i))*(S(i+1) - S(i-1)))...
    %    + r*S(i)/(S(i+1) - S(i));  % same as beta_central 
    %%backward alpha and beta formula
    alpha_backward(i) = alpha_parameter^2*S(i) / ((S(i) - S(i-1))*(S(i+1) - S(i-1))) ...
        - r * S(i) / (S(i+1) - S(i));
    beta_backward(i) = alpha_parameter^2*S(i)/((S(i+1) - S(i))*(S(i+1) - S(i-1))); 
    %alpha_backward(i) = sigma^2*S(i).^2 / ((S(i) - S(i-1))*(S(i+1) - S(i-1))) ...
    %    - r * S(i) / (S(i+1) - S(i));
    %beta_backward(i) = sigma^2*S(i).^2/((S(i+1) - S(i))*(S(i+1) - S(i-1)));
end

%% choosing parameter 

for i = 2:m-1
   if(alpha_central(i) >=0 && beta_central(i) >=0)
       alpha(i) = alpha_central(i);
       beta(i) = beta_central(i);
   elseif (alpha_forward(i) >=0 && beta_forward(i) >=0)
       alpha(i) = alpha_forward(i);
       beta(i) = beta_forward(i);
   else
       alpha(i) = alpha_backward(i);
       beta(i) = beta_backward(i);
   end
end


%% time step initialization
delt_sum = delt;
delt_old = delt;

%% CN-Rannacher time stepping

V_3 = zeros(m, N+1);
V_init_3 = max(K - S.^2, S.^2 - K)';
%V_init_3 = max(K-S, S-K)';
V_3(:,1) = V_init_3;
V_old_3 = V_init_3;
V_new_3 = V_old_3;
V_new3n = [];

for i = 1:2
    %% construct the m matrix
    %mid = delt_old.*(alpha + beta + r);
    %mid(m,:) = 0;
    %alpha(1,:) = [];
    %alpha = [alpha; 0];
    %beta(m,:) = [];
    %beta = [0; beta];

    M_matrix = [delt_old.*-alpha, delt_old.*(alpha + beta + r), delt_old.*-beta];
    M = spdiags(M_matrix, [-1,0,1], m-1, m);
    M = full([M;zeros(1,m)]);
    I = eye(m);
    t = 1;
    while t > tolerance
        pv = Large *(V_new_3 < Vs);
        PV = diag(pv);
        rhs1 = V_old_3 + PV * Vs; % RHS of the equation
        %AP = sparse(I + M + PV);
        AP = sparse(spdiags(ones(m,1),0,m,m)+M+PV);
        [L3,U3,P3,Q3] = lu(AP);
        V_new3n = Q3 * ((L3*U3)\(P3* rhs1));% compute (Vn+1)(k+1)
        t = max(abs(V_new3n - V_new_3)./(max(ones(m,1), abs(V_new3n))));
        V_new_3 = V_new3n;
    end
    %% implement timestep selector
    %bottom = [ones(m,1), abs(V_old_3), abs(V_new_3)];
    %top = abs(V_new_3 - V_old_3);
    %MaxRelChange = max(top./max(bottom,[],2));
    MaxRelChange = max(abs(V_new_3 -V_old_3)./(max(max(1,abs(V_new_3)),abs(V_old_3))));
    delt_new = (dnorm/MaxRelChange)*delt_old;
    delt_sum = delt_sum + delt_new;
    %sprintf("the value of delt_sum is %1d and the price value is : %1d",delt_sum, V_new_3(S == S0))
    
    delt_old = delt_new;
    
    
    V_old_3 = V_new_3;
end       
iteration = 2;

while delt_sum < T
    %% construct the m matrix
    %mid = delt_old.*(alpha + beta + r);
    %mid(m,:) = 0;
    %alpha(1,:) = [];
    %alpha = [alpha; 0];
    %beta(m,:) = [];
    %beta = [0; beta];

    M_matrix = [delt_old.*-alpha/2, delt_old.*(alpha + beta + r)/2, delt_old.*-beta/2];
    %M = full(spdiags(M_matrix, [-1,0,1], m, m));
    M = spdiags(M_matrix, [-1,0,1], m-1, m);
    M = full([M;zeros(1,m)]);
    I = eye(m);
    t =1;
    while t > tolerance
        pv = Large *(V_new_3 < Vs);
        PV = diag(pv);
        rhs2 = (I - M)* V_old_3 + PV*Vs; % RHS of the equation
        AP2 = sparse(spdiags(ones(m,1),0,m,m)+M+PV);
        %AP2 = sparse(I + 0.5*M + PV);
        [L4,U4,P4,Q4] = lu(AP2);
        V_new3n = Q4 * ((L4*U4)\(P4* rhs2));% compute (Vn+1)(k+1)
        t = max(abs(V_new3n - V_new_3)./(max(1, abs(V_new3n))));
         % compute relative change
         V_new_3 = V_new3n;

    end
    %bottom = [ones(m,1), abs(V_old_3), abs(V_new_3)];
    %top = abs(V_new_3 - V_old_3);
    %MaxRelChange = max(top./max(bottom,[],2));
    MaxRelChange = max(abs(V_new_3 -V_old_3)./(max(max(1,abs(V_new_3)),abs(V_old_3))));
    delt_new = (dnorm/MaxRelChange)*delt_old;
    delt_sum = delt_sum + delt_new;
    delt_old = delt_new;
    i = i + 1;
    %sprintf("the value of delt_sum is %d  and the price value is : %d ",delt_sum, V_new_3(S == S0))
    
    V_old_3 = V_new_3;
    iteration = iteration + 1;
end


delt_old = T-(delt_sum - delt_new);  %% last time step
if delt_old >0
    %mid = delt_old.*(alpha + beta + r);
    %mid(m,:) = 0;
    %alpha(1,:) = [];
    %alpha = [alpha; 0];
    %beta(m,:) = [];
    %beta = [0; beta];

    M_matrix = [delt_old.*-alpha/2, delt_old.*(alpha + beta + r)/2, delt_old.*-beta/2];
    %M = full(spdiags(M_matrix, [-1,0,1], m, m));
    M = spdiags(M_matrix, [-1,0,1], m-1, m);
    M = full([M;zeros(1,m)]);
    I = eye(m);
    t = 1
    while t > tolerance
  pv = Large *(V_new_3 < Vs);
        PV = diag(pv);
        rhs2 = (I - M)* V_old_3 + PV*Vs; % RHS of the equation
        AP2 = sparse(spdiags(ones(m,1),0,m,m)+M+PV);
        %AP2 = sparse(I + 0.5*M + PV);
        [L4,U4,P4,Q4] = lu(AP2);
        V_new3n = Q4 * ((L4*U4)\(P4* rhs2));% compute (Vn+1)(k+1)
        t = max(abs(V_new3n - V_new_3)./(max(1, abs(V_new3n))));
         % compute relative change
         V_new_3 = V_new3n;
    end    
    V_old_3 = V_new_3;
end
iteration = iteration + 1;
X1 = sprintf('The option value for fully implict in N = %d, S= %d, is: %d',N,length(S),V_new_3(S == S0));
disp(X1)
result = V_new_3(S == S0);
end
