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

figure(2);
S_CN_R_2 = S_CN_R(2:end);
plot(S_CN_R_2,delta)
xlabel('stock price')
ylabel('delta')

figure(3);
S_CN_R_3 = S_CN_R(3:end);
plot(S_CN_R_3,gamma)
xlabel('stock price')
ylabel('gamma')
