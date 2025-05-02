function P = ZCB(r,theta_r,g_r,c_r,t,T,U)
%Price at time t of a zero coupon bond with maturity T
Nsim = length(U);
psi = -(1-exp(-theta_r*(U*(T-t))))/theta_r;
P = exp(... 
    -r'*(1-exp(-theta_r*(T-t)))/theta_r...
    -g_r*sum(log(1-psi/c_r))*(T-t)/Nsim );
end