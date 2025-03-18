function X = simulateX(x0, alpha, sigma, t0, tEval, Nsim)
%
% INPUT
%
% x0 : initial condition
% alpha, sigma, : model parameters for the Gaussian OU
% t0: initial date of processes and seasonality
% tEval: Evaluation date
% Nsim : number of simulation
%
% OUTPUT
%
% X : simulation of the Gaussian OU process
    
    rng(0);
    X = x0*exp(-alpha*yearfrac(t0,tEval,3)) + sigma*exp(-alpha*yearfrac(t0,tEval,3))*sqrt( (exp(2*alpha*yearfrac(t0,tEval,3)) - 1) /(2*alpha) )*randn(Nsim,1);
end