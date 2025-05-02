function [P] = CDXOMC(t,T0,r,lambda,T,kappa,delta,U,Gr,Gl,...
    theta_r,theta_l,rho,c_r,c_l,c_t,g_r,g_l,g_t,N,L)
%r(T_0) and lambda(T_0) simulated via MonteCarlo, and CDX payoff computed
%accordingly.
dis = 0;
Nsim = length(Gr);
dt = (T0-t)/Nsim;
er = exp(-theta_r*((T0-t)-cumsum(ones(Nsim,1))*dt));
el = exp(-theta_l*((T0-t)-cumsum(ones(Nsim,1))*dt));
eY_r = (1-exp(-theta_r*((T0-t)-cumsum(ones(Nsim,1))*dt)))/theta_r;
eY_l = (1-exp(-theta_l*((T0-t)-cumsum(ones(Nsim,1))*dt)))/theta_l;
rT0 = r*exp(-theta_r*(T0-t))+Gr*er;
lT0 = lambda*exp(-theta_l*(T0-t))+(rho*Gr+Gl)*el;
Y_rT0 = r*(1-exp(-theta_r*(T0-t)))/theta_r+Gr*eY_r;
Y_lT0 = zeros(Nsim,L+1);%lambda*(1-exp(-theta_l*(T0-t)))/theta_l+(rho*Gr+Gl)*eY_l;
PP = zeros((N+1)*(L+1),Nsim);
check=0;
for s=1:Nsim
PP(:,s) = reshape((exp(-Y_rT0(s,:)-Y_lT0(s,:)'))',(N+1)*(L+1),1)...
         .*max(CDXPIDE(T0,T0,rT0(s,:),lT0(s,:),T,kappa,delta,U,...
         theta_r,theta_l,rho,c_r,c_l,c_t,g_r,g_l,g_t,N,L),0);
     if PP(:,s)>0
         check = check+1;
     end
end
P = sum(PP,2)/Nsim;
if dis==1
    [rT0,lT0]
    check
end
end

