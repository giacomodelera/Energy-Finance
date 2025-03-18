function prices = callPricesVG(alpha, sigma, beta, eta, A, B, C, nu, meanVG, sigmaVG, t0, tFinal, tau1,tau2, Nsim, step, x0, y0, discount, strike, meanPrice, stdPrice)
%
% INPUT
%
% alpha, sigma, beta, eta, A, B, C, nu, meanVG, sigmaVG : model parameters
% t0: initial date of processes and seasonality
% tFinal: Maturity date of the call
% tau1, tau2 : first and last monitoring dates
% step : discretization step in the Variance Gamma process
% Nsim : number of simulation
% x0, y0 : initial condition
% discount : discount factor
% strike: strike price(s)
% meanPrice, stdPrice : mean and stardard deviation of market prices
%
% OUTPUT
%
% prices : vector of call prices (1 x length(strike))

rng(0);

% swap prices
swaps = computeSwapsVG(alpha, sigma, beta, eta, A ,B, C, nu, meanVG, sigmaVG, t0, tFinal, tau1, tau2, step, Nsim, x0, y0);
swaps = meanPrice + stdPrice*swaps;

% Call prices
prices = discount*mean(max(swaps - strike,0));

end