function zero = eulereq(k_t,k_t1,k_t2,z_t,eps_z,rho,sigma)
% this function represent the stochastic euler equation that drives the
% capital accumultation (it holds in expected value)
% k_t, k_t1, k_t2 capital levels in period t,t+1, t+2
% z_t productivity, eps_z shock of the productivity
% rho and sigma parameter of the AR proces of the productivity

global alpha theta delta gamma beta

% the productivity follows a AR(1) process in log
z_t1=z_t^rho*exp(eps_z*sigma); 

% labor using the formula derived fro mthe F.O.C. w.r.t to labor (h_t)
h_t=((1-alpha)*z_t*k_t^alpha)^(1/(theta+alpha)); % h_t labor in period t
h_t1=((1-alpha)*z_t1*k_t1^alpha)^(1/(theta+alpha)); % h_t1 labor in pweiod t+1

% consumption in perio t and t+1
c_t=z_t * k_t^alpha * h_t^(1-alpha) + (1-delta)*k_t - k_t1; % c_t
c_t1=z_t1 * k_t1^alpha * h_t1^(1-alpha) + (1-delta)*k_t1 - k_t2; % c_t1
c_t=max(c_t,0.00001); % c_t cannot be negative
c_t1=max(c_t1,0.00001); % c_t1 cannot be negative

% marginal utility 
Up_t=(c_t - h_t^(1+theta) / (1+theta))^(-gamma); % marginal utility in period t
Up_t1=(c_t1 - h_t1^(1+theta) / (1+theta))^(-gamma); % marginal utility in period t+1

% decision rule
zero=Up_t-beta*Up_t1*(z_t1*alpha* (h_t/k_t)^(1-alpha) +1-delta);

end

