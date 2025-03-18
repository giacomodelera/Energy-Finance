clear all
clc
close all

%% CALIBRATION WITH ONLY THE GAUSSIAN OU

% same procedure of point 5 without Y process --> all parameters of the Y
% process equal to 0

% first sheet
data1 = readtable('DATA_DEEEX.xlsx', 'Sheet', 'Prices');
expiry_date = data1(:,1).Variables;
months = data1(:,2).Variables;
last = data1(:,3).Variables;
monitoringLength = ['M' 'M' 'M' 'M' 'M' 'M' 'M' 'Q' 'Q' 'Q' 'Q']';

% alpha : params(1)
% sigma : params(2)
% beta : 0
% eta : 0
% A : params (3)
% B : params (4)
% C : params (5)
% nu : 0
% meanVG : 0
% sigmaVG : 0

t0 = datetime(2024,1,1); % start date of processes and seasonality
today = datetime(2024,11,4);
step = 100;
Nsim = 500;
x0 = 1; y0 = 0;

tau1 = datetime.empty; % first monitoring date for each future
tau2 = datetime.empty; % last monitoring date for each future
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

%% CALIBRATION WITH STANDARDIZED PRICES

% alpha : params(1)
% sigma : params(2)
% beta : 0
% eta : 0
% A : params (3)
% B : params (4)
% C : params (5)
% nu : 0
% meanVG : 0
% sigmaVG : 0

err = @(params) (last-mean(last))/std(last) - gaussianOU_calibration(params(1), params(2), params(3), params(4), params(5), t0, today, tau1, tau2, Nsim, x0);
g = @(params) mean(err(params).^2);

%initial_guess = [alpha, sigma, A, B, C];
initial_guess = [0 0 0 0 0];

lb = [1e-8, 1e-8, -inf, -inf, -inf]; % Lower bounds
ub = [inf inf inf inf inf ]; % Upper bounds

fprintf("\n\n ----- CALIBRATION WITH STANDARDIZED PRICES ---- \n\n")
fprintf("OPTIMIZATION STARTED ...\n" )

tic
options = optimoptions('fmincon', 'Display','none','StepTolerance',1e-8, 'MaxIterations', 200);
global_coeffs_std = fmincon(g, initial_guess,[],[],[],[],lb,ub,[],options);
toc

fprintf("\n ... OPTIMIZATION FINISHED ! \n\n")

calibrated_price = mean(last)+(std(last)*gaussianOU_calibration(global_coeffs_std(1), global_coeffs_std(2), global_coeffs_std(3), global_coeffs_std(4), global_coeffs_std(5), t0, today, tau1, tau2, Nsim, x0));

figure()
plot(last, 'b','LineWidth',2)
hold on
grid on
plot(calibrated_price,'r', 'LineWidth',2)
title("GLOBAL CALIBRATION")
legend("Market prices", "Model prices", 'Location','southeast')

fprintf("RMSE : %d \n", sqrt(mean((last - calibrated_price).^2)))
fprintf("MAE : %d \n", mean(abs(last - calibrated_price)))

%% CALIBRATION WITH NON-STANDARDIZED PRICES

err = @(params) last - gaussianOU_calibration(params(1), params(2), params(3), params(4), params(5), t0, today, tau1, tau2, Nsim, x0);
g = @(params) mean(err(params).^2);

%initial_guess = [alpha, sigma, A, B, C];
initial_guess = [0 0 0 0 0];

lb = [0, 0, -inf, -inf, -inf]; % Lower bounds
ub = [inf inf inf inf inf ]; % Upper bounds

fprintf("\n\n ----- CALIBRATION WITH STANDARDIZED PRICES ---- \n\n")
fprintf("OPTIMIZATION STARTED ...\n" )

tic
options = optimoptions('fmincon', 'Display','none','StepTolerance',1e-8, 'MaxIterations', 200);
global_coeffs = fmincon(g, initial_guess,[],[],[],[],lb,ub,[],options);
toc

fprintf("\n ... OPTIMIZATION FINISHED ! \n\n")

calibrated_price = gaussianOU_calibration(global_coeffs(1), global_coeffs(2), global_coeffs(3), global_coeffs(4), global_coeffs(5), t0, today, tau1, tau2, Nsim, x0);

figure()
plot(last, 'b','LineWidth',2)
hold on
grid on
plot(calibrated_price,'r', 'LineWidth',2)
title("GLOBAL CALIBRATION")
legend("Market prices", "Model prices", 'Location','southeast')

fprintf("RMSE : %d \n", sqrt(mean((last - calibrated_price).^2)))
fprintf("MAE : %d \n", mean(abs(last - calibrated_price)))


%% CALL PRICING

tFinal = today + calmonths(6); % maturity date
while(isweekend(tFinal))
    tFinal = tFinal + caldays(1);
end

% Parameters
Nsim = 1e5;
alpha = global_coeffs_std(1); 
sigma = global_coeffs_std(2); 
A =  global_coeffs_std(3);
B =  global_coeffs_std(4);
C =  global_coeffs_std(5);

% 4Q25 swap --> tau1, tau2 and monitoring dates
tau1 = datetime(2025,10,1);     % tau1
tau2 = tau1 + calmonths(3) - caldays(1);    % tau2
monitoringDates = getMonitoringDays(tau1,tau2);

% discount factor @tFinal
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
calls = callPricesGaussian(alpha,sigma,A,B,C,t0,tFinal,tau1,tau2,Nsim,x0,disc_tFinal,strike, mean(last), std(last));



%% IMPLIED VOLATILITIES

futurePrice = last(end);
impl_vol = blsimpv(futurePrice, strike, r_6may,yearfrac(today,tFinal,3),calls);

figure()
plot(strike,impl_vol, '-o','LineWidth',2)
grid on
title("IMPLIED VOLATILITY")
 
