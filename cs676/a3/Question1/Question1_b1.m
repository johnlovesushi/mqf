alpha_parameter =0.8;    % constant related to alpha
r = .02;        % risk free rate
T = 1;          % time to expiry
K = 10;         % strike price
S0 = 10;        % initial stock price
N = 400;
delt = T/N;    % time interval
S = S4;

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

%% construct the m matrix

M_matrix = [delt.*-alpha,delt.*(alpha + beta + r) , delt.*-beta];
M = spdiags(M_matrix, [-1,0,1], m-1, m);
M = full([M;zeros(1,m)]);
I = eye(m);

A1 = sparse(I + M);     % create space matrix
[L1,U1,P1,Q1] = lu(A1); % lu factorization

M_matrix_CN = [delt.*-alpha/2, delt.*(alpha + beta + r)/2, delt.*-beta/2];
M_CN = spdiags(M_matrix_CN, [-1,0,1], m-1, m);
M_CN = full([M_CN;zeros(1,m)]);

A2 = sparse(spdiags(ones(m,1),0,m,m)+M_CN);     % create space matrix
[L2,U2,P2,Q2] = lu(A2); % lu factorization

V = zeros(m, N+1);
V_init = max(K - S.^2, S.^2 - K)';
V(:,1) = V_init;
V_old = V_init;
V_new = zeros(m,1);

for i = 1: N
            V_new = Q1 * ((L1*U1) \ P1) * V_old;
        V_old = V_new;
        V(:,i+1) = V_new; 
end
X1 = sprintf('The option value for fully implict in N = %d, S= %d, is: %s',N,length(S),V_new(S == S0));
disp(X1)
result1 = V_new(S == S0);

%% for CN time stepping 

A2 = sparse(I + 0.5*M);     % create space matrix
[L2,U2,P2,Q2] = lu(A2); % lu factorization

V_2 = zeros(m, N+1);
V_init_2 = max(K - S.^2, S.^2 - K)';
V_2(:,1) = V_init_2;
V_old_2 = V_init_2;
V_new_2 = zeros(m,1);

V_init_2 = max(K - S.^2, S.^2 - K)';
%V_init_2 = max(K - S, S - K)';
V_2(:,1) = V_init_2;
V_old_2 = V_init_2;
V_new_2 = zeros(m,1);

for i = 1: N
        Vt = (I - M_CN)*V_old_2;
        V_new_2 = Q2 * ((L2*U2) \ (P2*Vt)) ;
        V_old_2 = V_new_2;
        V_2(:,i+1) = V_new_2; 
end
X2 = sprintf('The option value for NC in N = %d, S= %d, is: %s',N,length(S),V_new_2(S == S0));
disp(X2)
result2 = V_new_2(S == S0);
%% CN-Rannacher time stepping


%M = spdiags(M_matrix, [-1,0,1], m-1, m);
%M = full([M;zeros(1,m)]);
V_3 = zeros(m, N+1);
V_init_3 = max(K - S.^2, S.^2 - K)';
%V_init_3 = max(K-S, S-K)'
V_3(:,1) = V_init_3;
V_old_3 = V_init_3;
V_new_3 = zeros(m,1);

%first two use fully implicit method
for i = 1: 2
        V_new_3 = Q1 * ((L1*U1) \ (P1 * V_old_3));
        V_old_3 = V_new_3;
        V_3(:,i+1) = V_new_3; 
end

for i = 3: N
        Vt_R = (I - M_CN)*V_old_3;
        V_new_3 = Q2 * ((L2*U2) \ (P2 * Vt_R));
        V_old_3 = V_new_3;
        V_3(:,i+1) = V_new_3; 
end
X3 = sprintf('The option value for NC-R in N = %d, S= %d, is: %s',N,length(S),V_new_3(S == S0));
disp(X3)
result3 = V_new_3(S == S0);