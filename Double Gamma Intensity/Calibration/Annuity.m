%% Annuity
% This function computes the annuity value ofthe index CDX swap at time t

function [A] = Annuity(t,T0,r,lambda,T,U,...
    theta_r,theta_l,rho,c_r,c_l,g_r,g_l,N,L)

Nsim = length(U);

M = length(T);
hvec = zeros((N+1)*(L+1),M);%Functions h_{\ell}

b3 = -rr(theta_r,T(1),T0);
b4 = -llambda(theta_l,T(1),T0);
xir_h = xi_r(T0-t,r,-1,b3,theta_r);
xil_h = xi_l(T0-t,lambda,0,b4,theta_l);
xi_h = reshape((xir_h+xil_h')',(N+1)*(L+1),1);
ppsi_r = psir(T0,t+U*(T0-t),-1,0,b3,b4,theta_r,theta_l,rho);
ppsi_l = psil(T0,t+U*(T0-t),0,b4,theta_l);
psi_r = psir(T(1),T0+U*(T(1)-T0),-1,-1,0,0,theta_r,theta_l,rho);
psi_l = psil(T(1),T0+U*(T(1)-T0),-1,0,theta_l);
hvec(:,1) = exp(  -g_r*sum(log(1-ppsi_r/c_r))*(T0-t)/Nsim+...
                  -g_l*sum(log(1-ppsi_l/c_l))*(T0-t)/Nsim+...
                  -g_r*sum(log(1-psi_r/c_r))*(T(1)-T0)/Nsim+...
                  -g_l*sum(log(1-psi_l/c_l))*(T(1)-T0)/Nsim+...
                   xi_h  );
for ell=2:M
    
    b3 = -rr(theta_r,T(ell),T0);
    b4 = -llambda(theta_l,T(ell),T0);
    xir_h = xi_r(T0-t,r,-1,b3,theta_r);
    xil_h = xi_l(T0-t,lambda,0,b4,theta_l);
    xi_h = reshape((xir_h+xil_h')',(N+1)*(L+1),1);
    ppsi_r = psir(T0,t+U*(T0-t),-1,0,b3,b4,theta_r,theta_l,rho);
    ppsi_l = psil(T0,t+U*(T0-t),0,b4,theta_l);
    psi_r = psir(T(ell),T0+U*(T(ell)-T0),-1,-1,0,0,theta_r,theta_l,rho);
    psi_l = psil(T(ell),T0+U*(T(ell)-T0),-1,0,theta_l);
    hvec(:,ell) = exp(  -g_r*sum(log(1-ppsi_r/c_r))*(T0-t)/Nsim+...
                        -g_l*sum(log(1-ppsi_l/c_l))*(T0-t)/Nsim+...
                        -g_r*sum(log(1-psi_r/c_r))*(T(ell)-T0)/Nsim+...
                        -g_l*sum(log(1-psi_l/c_l))*(T(ell)-T0)/Nsim+...
                        xi_h  );
end
TT = T-[T0,T(1:end-1)];
A = hvec*TT';
end

%% Other functions
function [xir] = xi_r(t,r,alpha1,alpha3,theta_r) % function \xi_r
xir=alpha1*r*rr(theta_r,t,0)+alpha3*r*exp(-theta_r*t);
end
function [xil] = xi_l(t,lambda,alpha2,alpha4,theta_l) % function \xi_l
xil=alpha2*lambda*llambda(theta_l,t,0)+alpha4*lambda*exp(-theta_l*t);
end
function [psiru] = psir(t,u,alpha1,alpha2,alpha3,alpha4,theta_r,theta_l,rho)
psiru = alpha1*rr(theta_r,t,u)+alpha2*rho*llambda(theta_l,t,u)...
    +alpha3*exp(-theta_r*(t-u))+alpha4*rho*exp(-theta_l*(t-u));
end
function [psilu] = psil(t,u,alpha2,alpha4,theta_l)
psilu = alpha2*llambda(theta_l,t,u)+alpha4*exp(-theta_l*(t-u));
end
function [rru] = rr(theta_r,t,u)
rru = (1-exp(-theta_r*(t-u)))/theta_r;
end
function [llambdau] = llambda(theta_lambda,t,u)
llambdau = (1-exp(-theta_lambda*(t-u)))/theta_lambda;
end