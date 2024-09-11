
%sigma = 0.15;
alpha_parameter =0.8;    % constant related to alpha
r = .02;                 % risk free rate
T = 1;                   % time to expiry
K = 10;                  % strike price
S0 = 10;                 % initial stock price
%N = 25;
delt = T/N;              % time interval
S = S4;
Large = 1e6;             % penalty coefficient
tolerance = 1/Large;     % toleranceeranceerance term
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

%% construct the m matrix
mid = delt.*(alpha + beta + r);
mid(m,:) = 0;
alpha(1,:) = [];
alpha = [alpha; 0];
beta(m,:) = [];
beta = [0; beta];

M_matrix = [delt.*-alpha, mid, delt.*-beta];
M = full(spdiags(M_matrix, [-1,0,1], m, m));
I = eye(m);


%% CN-Rannacher time stepping

V_3 = zeros(m, N+1);
V_init_3 = max(K - S.^2, S.^2 - K)';
%V_init_3 = max(K-S, S-K)';
V_3(:,1) = V_init_3;
V_old_3 = V_init_3;
V_new_3 = V_old_3;
V_new3n = [];

for i = 1:2
    for j = 1:10000
        pv = Large *(V_new_3 < Vs);
        PV = diag(pv);
        rhs1 = V_old_3 + PV * Vs; % RHS of the equation
        AP = sparse(I + M + PV);
        [L3,U3,P3,Q3] = lu(AP);
        V_new3n = Q3 * ((L3*U3)\P3)* rhs1;% compute (Vn+1)(k+1)
        t = max(abs(V_new3n - V_new_3)./(max(ones(m,1), abs(V_new3n))));
        V_new_3 = V_new3n;
        if t < tolerance
            break;
        end
    end
    V_old_3 = V_new_3;
end       


for i = 3:N
    for j = 1:10000
        pv = Large *(V_new_3 < Vs);
        PV = diag(pv);
        rhs2 = (I - 0.5*M)* V_old_3 + PV*Vs; % RHS of the equation
        AP2 = sparse(I + 0.5*M + PV);
        [L4,U4,P4,Q4] = lu(AP2);
        V_new3n = Q4 * ((L4*U4)\P4)* rhs2;% compute (Vn+1)(k+1)
        t = max(abs(V_new3n - V_new_3)./(max(ones(m,1), abs(V_new3n))));
         % compute relative change
         V_new_3 = V_new3n;
         if t <tolerance
             break;
         end
    end
    V_old_3 = V_new_3;
end
X1 = sprintf('The option value for fully implict in N = %d, S= %d, is: %s',N,length(S),V_new_3(S == S0));
disp(X1)
result1 = V_new_3(S == S0);
S_CN_R = S(S>=5 & S <=15)';
V_CN_R = V_new_3(S>=5 & S <=15);
n = length(S_CN_R);

delta = diff(V_CN_R)./diff(S_CN_R);
gamma=(((V_CN_R(3:n) - V_CN_R(2:n-1)) ./ (S_CN_R(3:n) - S_CN_R(2:n-1))) ...
    -((V_CN_R(2:n-1) - V_CN_R(1:n-2)) ./ (S_CN_R(2:n-1) - S_CN_R(1:n-2)))) ...
    ./ ((S_CN_R(3:n) - S_CN_R(1:n-2))/2);

figure(1);
plot(S_CN_R,V_CN_R);
xlabel('stock price')
ylabel('option value')
title('option value vs stock price in const timestep CN-R')

figure(2);
S_CN_R_2 = S_CN_R(2:end);
plot(S_CN_R_2,delta)
title('delta vs stock price in const timestep CN-R')
xlabel('stock price')
ylabel('delta')


%figure(3);
%S_CN_R_3 = S_CN_R(3:end);
%plot(S_CN_R_3,gamma)
%xlabel('stock price')
ylabel('gamma')