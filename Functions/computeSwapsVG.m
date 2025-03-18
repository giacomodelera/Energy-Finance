function swaps = computeSwapsVG(alpha, sigma, beta, eta, A ,B, C, nu, meanVG, sigmaVG, t0, tEval, tau1, tau2, step, Nsim, x0, y0)

    rng(0);

    % SIMULATE THE STOCHASTIC PART @tEval
    X_tEval = simulateX(x0, alpha, sigma, t0, tEval, Nsim);
    Y_tEval = simulateY(y0, beta, eta, nu ,meanVG, sigmaVG, t0, tEval, step, Nsim);

    % Characteristic exponent, Theta function and seasonality
    charExp_VG = @(t,T,u) -1/nu * log(1 - 1i*meanVG*nu*u + 0.5*sigmaVG^2*nu*u^2);
    Theta = @(t,T) exp( ( charExp_VG( t,T,-1i*eta*exp(-beta*yearfrac(t,T,3))) ) + sigma^2*( 1 - exp(-2*alpha*yearfrac(t,T,3)) )/(4*alpha) );
    Lambda = @(T) A*sin(2*pi*yearfrac(t0,T,3)) + B + C*yearfrac(t0,T,3);
    
    monitoringDates = getMonitoringDays(tau1, tau2);
    forwards = zeros(length(X_tEval),length(monitoringDates));
    for m = 1:length(monitoringDates)
        forwards(:,m) = Lambda(monitoringDates(m)) * (Theta(tEval,monitoringDates(m))) * exp( X_tEval*exp( -alpha*yearfrac(tEval,monitoringDates(m),3) ) + Y_tEval*exp(-beta*yearfrac(tEval,monitoringDates(m),3) ) ) ;
    end

    wSwap = [0 yearfrac(monitoringDates(1:end-1), monitoringDates(2:end), 3)]';

    swaps = forwards*wSwap;

end