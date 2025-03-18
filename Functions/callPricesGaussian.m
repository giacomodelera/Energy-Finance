function prices = callPricesGaussian(alpha, sigma, A, B, C, t0, tFinal, tau1,tau2, Nsim, x0, discount, strike, meanPrice, stdPrice)
%
% INPUT
%
% alpha, sigma, A, B, C : model parameters
% t0: initial date of processes and seasonality
% tFinal: Maturity date of the call
% tau1, tau2 : first and last monitoring dates
% Nsim : number of simulation
% x0 : initial condition
% discount : discount factor
% strike: strike price(s)
% meanPrice, stdPrice : mean and stardard deviation of market prices
%
% OUTPUT
%
% prices : vector of call prices (1 x length(strike))

rng(0);

% swap prices
swaps = computeSwaps_Gaussian(alpha, sigma, A ,B, C, t0, tFinal, tau1, tau2, Nsim, x0);
swaps = meanPrice + stdPrice*swaps;

% Call prices
prices = discount*mean(max(swaps - strike,0));

end