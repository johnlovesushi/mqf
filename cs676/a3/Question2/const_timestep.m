function [result] = const_timestep(S,N)

alpha_parameter =0.8;    % constant related to alpha
r = .02;                 % risk free rate
T = 1;                   % time to expiry
K = 10;                  % strike price
S0 = 10;                 % initial stock price

delt = T/N;              % time interval

Large = 1e6;             % penalty coefficient
tolerance = 1/Large;     % toleranceeranceerance term


Vs = max(K-S.^2, S.^2 - K)';

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

    %%forward alpha and beta formula
    alpha_forward(i) = alpha_parameter^2*S(i)/((S(i) - S(i-1))*(S(i+1) - S(i-1)));
    beta_forward(i) = alpha_parameter^2*S(i)/((S(i+1) - S(i))*(S(i+1) - S(i-1)))...
        + r*S(i)/(S(i+1) - S(i));  % same as beta_central

    %%backward alpha and beta formula
    alpha_backward(i) = alpha_parameter^2*S(i) / ((S(i) - S(i-1))*(S(i+1) - S(i-1))) ...
        - r * S(i) / (S(i+1) - S(i));
    beta_backward(i) = alpha_parameter^2*S(i)/((S(i+1) - S(i))*(S(i+1) - S(i-1))); 

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




%% CN-Rannacher time stepping

V_3 = zeros(m, N+1);
V_init_3 = max(K - S.^2, S.^2 - K)';

V_3(:,1) = V_init_3;
V_old_3 = V_init_3;
V_new_3 = V_old_3;
V_new3n = [];

for i = 1:2
    M_matrix = [delt.*-alpha, delt.*(alpha + beta + r), delt.*-beta];
    M = spdiags(M_matrix, [-1,0,1], m-1, m);
    M = full([M;zeros(1,m)]);
    I = eye(m);
    t = 1;
    while t > tolerance
        pv = Large *(V_new_3 < Vs);
        PV = diag(pv);
        rhs1 = V_old_3 + PV * Vs; % RHS of the equation

        AP = sparse(spdiags(ones(m,1),0,m,m)+M+PV);
        [L3,U3,P3,Q3] = lu(AP);
        V_new3n = Q3 * ((L3*U3)\(P3* rhs1));% compute (Vn+1)(k+1)
        t = max(abs(V_new3n - V_new_3)./(max(ones(m,1), abs(V_new3n))));
        V_new_3 = V_new3n;
    end
    V_old_3 = V_new_3;
end       


for i = 3:N
    M_matrix = [delt.*-alpha/2, delt.*(alpha + beta + r)/2, delt.*-beta/2];

    M = spdiags(M_matrix, [-1,0,1], m-1, m);
    M = full([M;zeros(1,m)]);
    I = eye(m);
    t =1;
    while t > tolerance
        pv = Large *(V_new_3 < Vs);
        PV = diag(pv);
        rhs2 = (I - M)* V_old_3 + PV*Vs; % RHS of the equation
        AP2 = sparse(spdiags(ones(m,1),0,m,m)+M+PV);

        [L4,U4,P4,Q4] = lu(AP2);
        V_new3n = Q4 * ((L4*U4)\(P4* rhs2));% compute (Vn+1)(k+1)
        t = max(abs(V_new3n - V_new_3)./(max(1, abs(V_new3n))));
         % compute relative change
         V_new_3 = V_new3n;
    end
    V_old_3 = V_new_3;
end
X1 = sprintf('The option value for fully implict in N = %d, S= %d, is: %s',N,length(S),V_new_3(S == S0));
disp(X1)
result = V_new_3(S == S0);
end