function swaps = computeSwaps_Gaussian(alpha, sigma, A ,B, C, t0, tEval, tau1, tau2, Nsim, x0)
%
% INPUT
%
% alpha, sigma, A, B, C, : model parameters
% t0: initial date of processes and seasonality
% tEval: Evaluation date
% tau1, tau2 : first and last monitoring dates
% Nsim : number of simulation
% x0 : initial condition
%
% OUTPUT
%
% swaps : vector of simulated swap prices @tEval (Nsim x 1)

    rng(0);

    % SIMULATE THE STOCHASTIC PART @tEval
    X_tEval = simulateX(x0, alpha, sigma, t0, tEval, Nsim);

    % Characteristic exponent, Theta function and seasonality
    Theta = @(t,T) exp( sigma^2*( 1 - exp(-2*alpha*yearfrac(t,T,3)) )/(4*alpha) );
    Lambda = @(T) A*sin(2*pi*yearfrac(t0,T,3)) + B + C*yearfrac(t0,T,3);
    
    monitoringDates = getMonitoringDays(tau1, tau2);

    % Simulated forward prices at monitoringDates maturity
    forwards = zeros(length(X_tEval),length(monitoringDates));
    for m = 1:length(monitoringDates)
        forwards(:,m) = Lambda(monitoringDates(m)) * (Theta(tEval,monitoringDates(m))) * exp( X_tEval*exp( -alpha*yearfrac(tEval,monitoringDates(m),3) ) ) ;
    end
    
    % forward weights
    wSwap = [0 yearfrac(monitoringDates(1:end-1), monitoringDates(2:end), 3)]';

    % swap price as weighted sum of forwards
    swaps = forwards*wSwap;

end