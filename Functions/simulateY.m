function Y = simulateY(y0, beta, eta, nu ,meanVG, sigmaVG, t0, tEval, step, Nsim)
%
% INPUT
%
% y0 : initial condition
% beta, eta, nu, meanVG, sigmaVG : model parameters for the additive OU
% t0: initial date of processes and seasonality
% tEval: Evaluation date
% step : discretization step in the Variance Gamma process
% Nsim : number of simulation
%
% OUTPUT
%
% Y : simulation of the additive OU process
    
    rng(0);
    dt = yearfrac(t0,tEval,3)/step;
    shapeGamma = 1 / nu * dt;
    rateGamma = 1 / nu;
    scaleGamma = 1/rateGamma;

    % increments of the gamma process
    increments = gamrnd(shapeGamma,scaleGamma,Nsim, step);
    
    X_VG = zeros(Nsim, step +1);

    for i = 1:Nsim
        % Gamma process
        G_t = cumsum(increments(i,:));
        
        % Surbordinated brownian motion
        W_Gamma = sqrt([0, G_t]) .* randn(1, step +1);
        
        % VG process
        X_VG(i, :) = meanVG * [0, G_t] + sigmaVG * W_Gamma;
    end

    % computation of Y process
    t_i = linspace(t0,tEval,step+1)';
    dt_i = cumsum(yearfrac(t_i(1:end-1), t_i(2:end),3)) - yearfrac(t_i(1),t_i(2),3)/2; 
    dX_VG = X_VG(:,2:end) - X_VG(:,1:end-1); % Variance gamma increments
    Y = y0*exp(-beta*yearfrac(t0,tEval,3)) + eta*exp(-beta*yearfrac(t0,tEval,3))*dX_VG*exp(beta*dt_i)*(tEval>t0);
end