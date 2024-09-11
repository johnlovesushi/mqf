function [Var,CVaR]=dVaRCVaR(PL, beta)
    %% this function is used to calculate the VaR and CVaR value
    %% It reqiure an input of Profit and Loss and beta
    sorted_PL = sort(PL);       % sort P&L in increasing order
    index = floor((1-beta)*length(PL));    % find the index of M simulated P&L for beta (s.t. index/M <= 1-beta)
    Var = sorted_PL(index);    % compute VaR
    CVaR = mean(sorted_PL(1:index));   %compte CVaR
end