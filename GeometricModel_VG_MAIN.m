clear all
clc
close all

%% Import future prices

% first sheet
data1 = readtable('DATA_DEEEX.xlsx', 'Sheet', 'Prices');
expiry_date = data1(:,1).Variables;
months = data1(:,2).Variables;
last = data1(:,3).Variables;
monitoringLength = ['M' 'M' 'M' 'M' 'M' 'M' 'M' 'Q' 'Q' 'Q' 'Q']';

t0 = datetime(2024,1,1); % Start date of the processes and seasonality
today = datetime(2024,11,4);
step = 100; 
Nsim = 500;
x0 = 1; y0 = 1;

tau1 = datetime.empty; % first monitoring dates for each future
tau2 = datetime.empty; % last monitoring dates for each future
for n = 1:length(last)
    yy = year(expiry_date(n));
    mm = month(expiry_date(n));
    tau1(n) = datetime(yy,mm+1,1);
    if (monitoringLength(n) == 'M')
        tau2(n) = tau1(n) + calmonths(1) - caldays(1);
    else
        tau2(n) = tau1(n) + calmonths(3) - caldays(1);
    end
end

%% Calibrate 𝛼,𝜎,𝛽,𝜂 𝐴,𝐵,𝐶 and the Variance Gamma parameters

% alpha : params(1)
% sigma : params(2)
% beta : params(3)
% eta : params (4)
% A : params (5)
% B : params (6)
% C : params (7)
% nu : params (8) --> variance rate Gamma process
% meanVG : params(9)
% sigmaVG : params(10)

%% CALIBRATION WITH STANDARDIZED PRICES

err = @(params) (last-mean(last))/std(last) - VG_Calibration(params(1), params(2), params(3), params(4), params(5), params(6), params(7), params(8), params(9), params(10), t0, today, tau1, tau2, step, Nsim, x0, y0);
g = @(params) mean(err(params).^2); % minimizing the mse

%initial_guess = [alpha, sigma, beta, eta, A, B, C, nu, meanVG, sigmaVG];
initial_guess = [0 0 0 0 0 0 0 1 0 0];

lb = [1e-8, 1e-8, 1e-8, 1e-8, -inf, -inf, -inf, 1e-8, -inf, 1e-8]; % Lower bounds
ub = [inf inf inf inf inf inf inf inf inf inf]; % Upper bounds

fprintf("\n\n ----- CALIBRATION WITH STANDARDIZED PRICES ---- \n\n")

fprintf("OPTIMIZATION STARTED ...\n" )
tic
options = optimoptions('fmincon', 'Display','none','StepTolerance',1e-8);
global_coeffs_std = fmincon(g, initial_guess,[],[],[],[],lb,ub,[],options);
toc
fprintf("... OPTIMIZATION FINISHED ! \n\n")

calibrated_price = mean(last)+(std(last)*VG_Calibration(global_coeffs_std(1), global_coeffs_std(2), global_coeffs_std(3), global_coeffs_std(4), global_coeffs_std(5), global_coeffs_std(6), global_coeffs_std(7), global_coeffs_std(8), global_coeffs_std(9), global_coeffs_std(10), t0, today, tau1, tau2, step, Nsim, x0, y0));

figure("name","Calibration with Standardized Prices")
plot(last, 'b','LineWidth',2)
hold on
grid on
plot(calibrated_price,'r', 'LineWidth',2)
title("GLOBAL CALIBRATION")
legend("Market prices", "Model prices", 'Location','southeast')

fprintf("RMSE : %d \n", sqrt(mean((last - calibrated_price).^2)))
fprintf("MAE : %d \n", mean(abs(last - calibrated_price)))

%% CALIBRATION WITH NON-STANDARDIZED PRICES
err = @(params) last - VG_Calibration(params(1), params(2), params(3), params(4), params(5), params(6), params(7), params(8), params(9), params(10), t0, today, tau1, tau2, step, Nsim, x0, y0);
g = @(params) mean(err(params).^2);

%initial_guess = [alpha, sigma, beta, eta, A, B, C, nu, meanVG, sigmaVG];
initial_guess = [0 0 0 0 0 0 0 1 0 0];

lb = [1e-8, 1e-8, 1e-8, 1e-8, -inf, -inf, -inf, 1e-8, -inf, 1e-8]; % Lower bounds
ub = [inf inf inf inf inf inf inf inf inf inf]; % Upper bounds


fprintf("\n\n ----- CALIBRATION WITH NO STANDARDIZED PRICES ---- \n\n")

fprintf("OPTIMIZATION STARTED ...\n")
tic
options = optimoptions('fmincon', 'Display','none','StepTolerance',1e-8);
global_coeffs = fmincon(g, initial_guess,[],[],[],[],lb,ub,[],options);
toc
fprintf("... OPTIMIZATION FINISHED ! \n\n")


calibrated_price_ = VG_Calibration(global_coeffs(1), global_coeffs(2), global_coeffs(3), global_coeffs(4), global_coeffs(5), global_coeffs(6), global_coeffs(7), global_coeffs(8), global_coeffs(9), global_coeffs(10), t0, today, tau1, tau2, step, Nsim, x0, y0);

figure("name","Calibration with NO Standardized Prices")
plot(last, 'b' ,'LineWidth',2)
hold on
grid on
plot(calibrated_price_,'r', 'LineWidth',2)
title("GLOBAL CALIBRATION")
legend("Market prices", "Model prices", 'Location','southeast')

fprintf("MSE : %d \n", sqrt(g(global_coeffs)))
fprintf("MAE : %d \n", mean(abs(last - calibrated_price_)))

%% Since the calibration from standardized prices performed better, it was decided to use the coeffients 
%% from that calibration and then rescale all the results with the mean and the standard deviation of market prices

%% Data for call prices

tFinal = today + calmonths(6);  % maturity date
while(isweekend(tFinal))
    tFinal = tFinal + caldays(1);
end

Step = 100; Nsim = 1e4;
alpha = global_coeffs_std(1); 
sigma = global_coeffs_std(2); 
beta = global_coeffs_std(3); 
eta = global_coeffs_std(4); 
A =  global_coeffs_std(5);
B =  global_coeffs_std(6);
C =  global_coeffs_std(7);
nu=  global_coeffs_std(8);
meanVG =  global_coeffs_std(9);
sigmaVG =  global_coeffs_std(10);

% 4Q25 swap --> tau1, tau2 and monitoring dates
tau1 = datetime(2025,10,1);     % tau1
tau2 = tau1 + calmonths(3) - caldays(1);    % tau2
monitoringDates = getMonitoringDays(tau1,tau2);

% Discount factor @tFinal
data3 = readtable('DATA_DEEEX.xlsx', 'Sheet', 'discounts');
disc_6apr = data3{2,8};
t1 = datetime(data3{1,8}, 'ConvertFrom', 'excel');
r_6apr = - log(disc_6apr)/yearfrac(today,t1,3);

disc_8may = data3{2,9};
t2 = datetime(data3{1,9}, 'ConvertFrom', 'excel');
r_8may = - log(disc_8may)/yearfrac(today,t2,3);

r_6may = interp1(datenum([t1 t2]), [r_6apr r_8may], datenum(tFinal), 'linear');

disc_tFinal = exp(-r_6may*yearfrac(today,tFinal,3));

% Strikes
strike = 400:10:600;

% Call prices
calls = callPricesVG(alpha,sigma,beta,eta,A,B,C,nu,meanVG,sigmaVG,t0,tFinal,tau1,tau2,Nsim,Step,x0,y0,disc_tFinal,strike, mean(last), std(last));


figure()
plot(strike,calls, '-o','Color',"r");
grid on
title("CALL PRICES")



%% Implied volatilities

mktPrice = last(end);
impl_vol = blsimpv(mktPrice, strike, r_6may, yearfrac(today,tFinal,3), calls);

figure()
plot(strike,impl_vol, '-o','LineWidth',2)
grid on
title("IMPLIED VOLATILITY")
