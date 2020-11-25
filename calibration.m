function zero = calibration(parameters,k_star,k1,z1,eps_z,a,b,D2,D3)
% this function is used to calibrate the model to pick rho (persistance of 
% the productivity) and sigma (volatility of the shocks) so the the model
% replicates the volatility and the autocorrelation of log output mesured by 
% kydland and Prescott
% use FSOLVE to solve for rho and sigma
% k_t determistic s.s. capital
% z_t productivity
% eps_z vector of the shocks
% a, b coeffients of the log linear solutiom
% D2, D3 derivates of the delta function w.r.t k_t+1 and k_t+2
global alpha T delta


rho=parameters(1); % persistance of volatility in log
sigma=parameters(2); % std dev of the shocks of the volatilty
h=10^-7; % increment to compute numerically the derivative
D4=imag(eulereq(k_star,k_star,k_star,1+1i*h,0,rho,sigma))/h; % derivateive of the euler equation w.r.t z_t
f=-D4/(k_star*(D2+D3*(b+rho))); % coefficient of the log-linear solution

% creates matrices to store results
log_k=zeros(T,1);
log_z=zeros(T,1);
log_h=zeros(T,1);
log_output=zeros(T,1);

% assign starting points
log_k(1)=log(k1);
log_z(1)=log(z1);
log_h(1)=log(labor(k1,z1));
log_output(1)=log_z(1)+alpha*log_k(1)+(1-alpha)*log_h(1);

% compute the cycle
for t=2:T
    log_k(t)=max(a+b*log_k(t-1)+f*log_z(t-1),log(1-delta)+log_k(t-1)); % I>=0, k_t+1 cannot be lower than (1-delta)*k_t
    log_z(t)=rho*log_z(t-1)+sigma*eps_z(t); % Ar(1) process of volatility
    log_h(t)=log(labor(exp(log_k(t)),exp(log_z(t)))); % log labor
    log_output(t)=log_z(t)+alpha*log_k(t)+(1-alpha)*log_h(t); %log(z*k^alpha h^(1-alpha))
end

% compute volatility and standard deviation
volatility_y=std(log_output); % std dev od log output
autocorrelation_y=autocorr(log_output); % autocorrelation of log autput

% impose Kydland and Prescott figures
zero1=volatility_y-3.5;
zero2=autocorrelation_y(2)-0.66;

% output of the function
zero=[zero1; zero2];
end