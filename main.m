% main.m
% this script compute the business cycle using a log-linear solution
% the productivity follows an AR(1) process in log



clear all
clc


% parameterization
global alpha beta delta theta gamma T
alpha=0.3; % parameter of Cobb-Douglass production function (CRS)
beta=1.04^-1; % discounting factor
delta=0.08; % depreciation rate
theta=0.6; % parameter of the disutility of labor
gamma=1.5; % parameter of the utility function
T=100;
% fix the seed of the randon number generator
rng(1) 

% draw random numbers from the normal distribution to simulate productivity
% shocks
eps_z=normrnd(0,1,T,1);

% compute deterministic steady state
k_star=((1/beta-1+delta)^-1 *alpha*(1-alpha)^((1-alpha)/(theta+alpha)))^((theta+alpha)/(theta*(1-alpha))); % capital
h_star=labor(k_star,1); % labor (in the deterministic staedy state the productivity (z) is assumed = 1)
z_star=1; % productivity in the steady state = 1 in every period

% compute derivates of euler equation 
% (to compute D4 (derivative w.r.t. z_t)we need  the value of rho)---> D4 
% computed after the calibration
% to compute D1, D2 and D3, the values rho and sigma are not necessary 

% standard method (numerical differentiation)
h=10^-5; % increment
D1=(eulereq(k_star+h,k_star,k_star,z_star,0,1,1)-eulereq(k_star-h,k_star,k_star,z_star,0,1,1))/(2*h);
D2=(eulereq(k_star,k_star+h,k_star,z_star,0,1,1)-eulereq(k_star,k_star-h,k_star,z_star,0,1,1))/(2*h);
D3=(eulereq(k_star,k_star,k_star+h,z_star,0,1,1)-eulereq(k_star,k_star,k_star-h,z_star,0,1,1))/(2*h);


% complex step differentiation
h_c=10^-7; % increment (it can be smaller than h)
D1c=imag(eulereq(k_star+i*h_c,k_star,k_star,z_star,0,1,1))/h_c;
D2c=imag(eulereq(k_star,k_star+i*h_c,k_star,z_star,0,1,1))/h_c;
D3c=imag(eulereq(k_star,k_star,k_star+i*h_c,z_star,0,1,1))/h_c;


% compute the coefficient of the log-liner solution

% b
% create the vectors of the polynomial used to solve for b
p=[D3c D2c D1c];
b=roots(p);
b=min(b); % Pick the stable root (0<b<1), the root>1 is unstable
% a
a=(1-b)*log(k_star);


% calibration

% pick rho and sigma to replicates the values of volatility and
% autocorrelation of log output mesured by kydland and Prescott
% set as starting point the deterministic S.S (k=k_star; z_star=1)
% use derivates computed through the complex step differentiation

x0=[0.5; 1]; % guess for the fsolve
solution=fsolve(@(parameters) calibration(parameters,k_star,k_star,1,eps_z,a,b,D2c,D3c),x0);
rho=solution(1); % autocorrelation of productivity in logs
sigma=solution(2); % std deviation of the shocks of the productivity


% compute the derivate of the euler equation w.r.t. to z_t using the values
% of rho and sigma obtained with the calibration

% standard method (numerical differentiation)
D4=(eulereq(k_star,k_star,k_star,z_star+h,0,rho,sigma)-eulereq(k_star,k_star,k_star,z_star-h,0,rho,sigma))/(2*h);

% complex step differentiation
D4c=imag(eulereq(k_star,k_star,k_star,z_star+1i*h_c,0,rho,sigma))/h_c;

% compute coeffient f of the log-linear solution (using derivates computed 
% through the complex step differentiation)
f=-D4c * (k_star)^-1 / (D2c+D3c*(b+rho)) ;

% compute the cycle

% create matrices to store results in
log_k=zeros(T,1);
log_z=zeros(T,1);

% assign starting points (= to their deterministic s.s. levels)
log_k(1)=log(k_star);
log_z(1)=log(z_star);

% compute the paths
for t=2:T
    log_k(t)=a+b*log_k(t-1)+f*log_z(t-1);
    log_k(t)=max(log_k(t),log(1-delta)+log_k(t-1)); % Investment>=0 (k_t+1>=(1-delta)*k_t)
    log_z(t)=rho*log_z(t-1)+sigma*eps_z(t);
end

% take the antilog to every series
k_t=exp(log_k);
z_t=exp(log_z);

% compute labor and output
h_t=labor(k_t,z_t);
output_t=z_t.*k_t.^alpha .*h_t.^(1-alpha);

% take the log of the series of c and Y
log_h=log(h_t);
log_output=log(output_t);


% compute investment and consumption

% create matrix to store results
i_t=zeros(T,1);
c_t=zeros(T,1);

% compute the paths
for t=1:T-1
    i_t(t)=k_t(t+1)-(1-delta)*k_t(t);
    i_t(t)=max(i_t(t),0.00001); % I cannot be negative
    c_t(t)=output_t(t)-i_t(t); % C+I=Y-->C=Y-I
    c_t(t)=max(c_t(t),0.00001); % consumption cannot be negative
end

% assign final values assuming that the values in period T+1 are the same as in T
i_t(T)=k_t(T)-(1-delta)*k_t(T);
c_t(T)=output_t(T)-i_t(T);

% compute the log of the series of c and i
log_i=log(i_t);
log_c=log(c_t);

% compute statistics of the log series

% volatility (standard deviation)
volatility_y=std(log_output);
volatility_c=std(log_c);
volatility_i=std(log_i);
volatility_h=std(log_h);
volatility_z=std(log_z);

% correlation with output
corr_yvsy=1;
corr_yvsc=corrcoef(log_output,log_c);
corr_yvsc=corr_yvsc(2);
corr_yvsi=corrcoef(log_output,log_i);
corr_yvsi=corr_yvsi(2);
corr_yvsh=corrcoef(log_output,log_h);
corr_yvsh=corr_yvsh(2);
corr_yvsz=corrcoef(log_output,log_z);
corr_yvsz=corr_yvsz(2);

% autocorrelation 1st lag
autocorrelation_y=autocorr(log_output);
autocorrelation_y=autocorrelation_y(2);
autocorrelation_c=autocorr(log_c);
autocorrelation_c=autocorrelation_c(2);
autocorrelation_i=autocorr(log_i);
autocorrelation_i=autocorrelation_i(2);
autocorrelation_h=autocorr(log_h);
autocorrelation_h=autocorrelation_h(2);
autocorrelation_z=autocorr(log_z);
autocorrelation_z=autocorrelation_z(2);

% display statitics
display("The statistics of the model using log series are:")
display("Volatility (std dev)")
display("Output Consumption Investment Labor Productivity")
display([volatility_y volatility_c volatility_i volatility_h volatility_z])

display("Correlation with output")
display("Output Consumption Investment Labor Productivity")
display([corr_yvsy corr_yvsc corr_yvsi corr_yvsh corr_yvsz])

display("Autocorrelation (1st lag)")
display("Output Consumption Investment Labor Productivity")
display([autocorrelation_y autocorrelation_c autocorrelation_i autocorrelation_h autocorrelation_z])


% compute statistics of the series (without logs)

% volatility (standard deviation)
volatility_ynolog=std(output_t);
volatility_cnolog=std(c_t);
volatility_inolog=std(i_t);
volatility_hnolog=std(h_t);
volatility_znolog=std(z_t);

% correlation with output
corr_yvsynolog=1;
corr_yvscnolog=corrcoef(output_t,c_t);
corr_yvscnolog=corr_yvscnolog(2);
corr_yvsinolog=corrcoef(output_t,i_t);
corr_yvsinolog=corr_yvsinolog(2);
corr_yvshnolog=corrcoef(output_t,h_t);
corr_yvshnolog=corr_yvshnolog(2);
corr_yvsznolog=corrcoef(output_t,z_t);
corr_yvsznolog=corr_yvsznolog(2);

% autocorrelation 1st lag
autocorrelation_ynolog=autocorr(output_t);
autocorrelation_ynolog=autocorrelation_ynolog(2);
autocorrelation_cnolog=autocorr(c_t);
autocorrelation_cnolog=autocorrelation_cnolog(2);
autocorrelation_inolog=autocorr(i_t);
autocorrelation_inolog=autocorrelation_inolog(2);
autocorrelation_hnolog=autocorr(h_t);
autocorrelation_hnolog=autocorrelation_hnolog(2);
autocorrelation_znolog=autocorr(z_t);
autocorrelation_znolog=autocorrelation_znolog(2);

% display statitics
display("The statistics of the model are (no log):")
display("Volatility (std dev)")
display("Output Consumption Investment Labor Productivity")
display([volatility_ynolog volatility_cnolog volatility_inolog volatility_hnolog volatility_znolog])

display("Correlation with output")
display("Output Consumption Investment Labor Productivity")
display([corr_yvsynolog corr_yvscnolog corr_yvsinolog corr_yvshnolog corr_yvsznolog])

display("Autocorrelation (1st lag)")
display("Output Consumption Investment Labor Productivity")
display([autocorrelation_ynolog autocorrelation_cnolog autocorrelation_inolog autocorrelation_hnolog autocorrelation_znolog])


% plot Capital and Productivity

figure(1)
subplot(1,2,1)
plot(log_k)
xlabel("time")
ylabel("log Capital")
title("Capital (in log) vs time")
subplot(1,2,2)
plot(k_t)
xlabel("time")
ylabel("Capital")
title("Capital vs time")


figure(2)
subplot(1,2,1)
plot(log_z)
xlabel("time")
ylabel("log Productivity")
title("Productivity (in log) vs time")
subplot(1,2,2)
plot(log_z)
xlabel("time")
ylabel("log Productivity")
title("Productivity (in log) vs time")


% display mean of capital stock and Productivity
% in log
display("Mean of capital stock in log")
display(mean(log_k))
display("Mean of the Productivity level in log")
display(mean(log_z))

% no log
display("Mean of capital stock")
display(mean(k_t))
display("Mean of the Productivity level")
display(mean(z_t))


% display solution
display("the coefficient of the log-linear solution are:")
display(a)
display(b)
display(f)


display("the parameters of the AR process of the productivity are:")
display("persistance")
display(rho)
display("volatility")
display(sigma)

% display derivates of the eurler equation
display("Derivates computed trhough the standard method")
display([D1 D2 D3 D4])
display("Derivates computed trhough the complex step differentiation")
display([D1c D2c D3c D4c])


% chech if the two types match
% create 2 vectors with the derivates
stdmethod=[D1 D2 D3 D4];
complexstep=[D1c D2c D3c D4c];
if norm(stdmethod-complexstep)>10^-4
    display("ATTENTION: they don't match!")
else 
    display("the derivates match")
end


% plot log series of capital, productivity, hours worked and output


figure(3)
title("log series of capital, productivity, hours worked and output ")
subplot(2,2,1)
plot(log_k)
xlabel("Time")
ylabel("log k")
subplot(2,2,2)
plot(log_z)
xlabel("Time")
ylabel("log z")
subplot(2,2,3)
plot(log_h)
xlabel("Time")
ylabel("log h")
subplot(2,2,4)
plot(log_output)
xlabel("time")
ylabel("log Y")

figure(4)
plot(log_h())
hold on
plot(log_output)
plot(log_k)
title("Output, capital and hours worked (in log)")
xlabel("time")
legend("h","Y","k")
hold off

figure(5)
plot(log_c(1:50))
xlabel("time")
ylabel("log C")
title("Consumption(log) vs time")

figure(6)
plot(log_i)
hold on
plot(log_k)
hold off
title("Investment and Capital (in log)")
xlabel("time")

% HP filter
figure(7)
hpfilter(log_k,1000);


figure(8)
hpfilter(log_output,1000);