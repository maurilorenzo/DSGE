function h = labor(k,z)
% this function computes the value of labor given the values of capital (k) 
% and productivity (z)
% fomrula derived from the F.O.C. w.r.t h_t
global alpha theta
% h_t
h=((1-alpha)*z.*k.^alpha).^(1/(theta+alpha));
end

