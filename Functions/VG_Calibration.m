function swaps = VG_Calibration(alpha, sigma, beta, eta, A ,B, C, nu, meanVG, sigmaVG, t0, tEval, tau1, tau2, step, Nsim, x0, y0)
%
% INPUT
%
% alpha, sigma, beta, eta, A, B, C, nu, meanVG, sigmaVG : model parameters
% t0: initial date of processes and seasonality
% tEval: Evaluation date
% tau1, tau2 : first and last monitoring dates
% step : discretization step in the Variance Gamma process
% Nsim : number of simulation
% x0, y0 : initial condition
%
% OUTPUT
%
% swaps : vector 11x1 of swap prices 

    rng(0);

    % SIMULATE THE STOCHASTIC PART @tEval
    X_tEval = simulateX(x0, alpha, sigma, t0, tEval, Nsim);
    Y_tEval = simulateY(y0, beta, eta, nu ,meanVG, sigmaVG, t0, tEval, step, Nsim);

    % Characteristic exponent of the Variance Gamma, Theta function and seasonality
    charExp_VG = @(t,T,u) ((-1/nu) * log(1 - 1i*meanVG*nu*u + 0.5*sigmaVG^2*nu*u^2));
    Theta = @(t,T) exp( charExp_VG( t,T,-1i*eta*exp(-beta*yearfrac(t,T,3)) ) + sigma^2*( 1 - exp(-2*alpha*yearfrac(t,T,3)) )/(4*alpha) );
    Lambda = @(T) A*sin(2*pi*yearfrac(t0,T,3)) + B + C*yearfrac(t0,T,3);
    
    % swap prices
    index = 1;
    swaps = zeros(length(tau1),1);
    for n = 1 : length(tau1)

        %forward prices
        monitoringDates = getMonitoringDays(tau1(n), tau2(n));
        forwards = zeros(length(X_tEval),length(monitoringDates));
        for m = 1:length(monitoringDates)
            forwards(:,m) = Lambda(monitoringDates(m)) *  (Theta(tEval,monitoringDates(m))) * exp( X_tEval*exp( -alpha*yearfrac(tEval,monitoringDates(m),3) ) + Y_tEval*exp(-beta*yearfrac(tEval,monitoringDates(m),3) ) ) ;
        end

        % integral weights
        wSwap = [0 yearfrac(monitoringDates(1:end-1), monitoringDates(2:end), 3)]';

        % swap price as weighted sum of forwards
        swaps(index) = mean(forwards*wSwap);
        index = index + 1;
     end


end