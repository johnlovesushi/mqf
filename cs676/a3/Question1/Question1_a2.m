K = 10;
N = 25;
S = [0:0.1*K:0.4*K,...  %input S value
    0.45*K:0.05*K:0.8*K,...
    0.82*K:0.02*K:0.9*K,...
    0.91*K:0.01*K:1.1*K,...
    1.12*K:0.02*K:1.2*K,...
    1.25*K:.05*K:1.6*K,...
    1.7*K:0.1*K:2*K,...
    2.2*K, 2.4*K, 2.8*K,...
    3.6*K, 5*K, 7.5*K, 10*K];
% for node 62 and timestep 25
L = length(S);
[V,V_CN,V_CNR] = stepping(S,N);  

% for node 123 and timestep 50
N1 = 50;
S1 = movmean(S,2);
S1 = [S,S1];
S1 = sort(S1);
S1 = S1(2:end);
L1 = length(S1);
[V1,V_CN1,V_CNR1] = stepping(S1,N1);

% for node 245 and timestep 100
N2 = 100;
S2 = movmean(S1,2);
S2 = [S1,S2];
S2 = sort(S2);
S2 = S2(2:end);
L2 = length(S2);
[V2,V_CN2,V_CNR2] = stepping(S2,N2);

% for node 489 and timestep 200
N3 = 200;
S3 = movmean(S2,2);
S3 = [S2,S3];
S3 = sort(S3);
S3 = S3(2:end);
L3 = length(S3);
[V3,V_CN3,V_CNR3] = stepping(S3,N3);

% for node 489 and timestep 200
N4 = 400;
S4 = movmean(S3,2);
S4 = [S3,S4];
S4 = sort(S4);
S4 = S4(2:end);
L4 = length(S4);
[V4,V_CN4,V_CNR4] = stepping(S4,N4);

T_1 =table([N;N1;N2;N3;N4],...
           [L;L1;L2;L3;L4],...
           [V;V1;V2;V3;V4],...
           [NaN;V1-V;V2-V1;V3-V2;V4-V3],...
           [NaN;NaN;(V1-V)/(V2-V1);(V2-V1)/(V3-V2);(V3-V2)/(V4-V3)]);
T_1.Properties.VariableNames ={'Timesteps','Node','Value','Change','Ratio'};

T_2 =table([N;N1;N2;N3;N4],...
           [L;L1;L2;L3;L4],...
           [V_CN;V_CN1;V_CN2;V_CN3;V_CN4],...
           [NaN;V_CN1-V_CN;V_CN2-V_CN1;V_CN3-V_CN2;V_CN4-V_CN3],...
           [NaN;NaN;(V_CN1-V_CN)/(V_CN2-V_CN1);(V_CN2-V_CN1)/(V_CN3-V_CN2);(V_CN3-V_CN2)/(V_CN4-V_CN3)]);
T_2.Properties.VariableNames ={'Timesteps','Node','Value','Change','Ratio'};


T_3 =table([N;N1;N2;N3;N4],...
           [L;L1;L2;L3;L4],...
           [V_CNR;V_CNR1;V_CNR2;V_CNR3;V_CNR4],...
           [NaN;V_CNR1-V_CNR;V_CNR2-V_CNR1;V_CNR3-V_CNR2;V_CNR4-V_CNR3],...
           [NaN;NaN;(V_CNR1-V_CNR)/(V_CNR2-V_CNR1);(V_CNR2-V_CNR1)/(V_CNR3-V_CNR2);(V_CNR3-V_CNR2)/(V_CNR4-V_CNR3)]);
T_3.Properties.VariableNames ={'Timesteps','Node','Value','Change','Ratio'};