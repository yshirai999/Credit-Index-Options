%% Forward CDX Price
% N.B: 
% (i) Allows defaults before maturity of forward,
% (ii) Excludes FEP (whose expected value can be computed using the EFEP
%      function and can be added to the price obtained here).
function [PI] = CDX(t,T0,r,lambda,T,kappa,delta,U,...
    theta_r,theta_l,rho,c_r,c_l,c_t,g_r,g_l,g_t,N,L)

Nsim = length(U);

M = length(T);
gvec = zeros((N+1)*(L+1),M); %Functions g_{\ell}
hvec = zeros((N+1)*(L+1),M); %Functions h_{\ell}
Id_M = ones(M,1);

alpha3 = -rr(theta_r,T(1),T0);
xir_g = xi_r(T0-t,r,-1,alpha3,theta_r);
xil_g = xi_l(T0-t,lambda,-1,0,theta_l);
xi_g = reshape((xir_g+xil_g')',(N+1)*(L+1),1);
psi_r = psir(T0,t+U*(T0-t),-1,-1,alpha3,0,theta_r,theta_l,rho);
psi_l = psil(T0,t+U*(T0-t),-1,0,theta_l);
P = exp(  -g_r*sum(log(1+rr(theta_r,T(1),T0+U*(T(1)-T0))/c_r))*...
           (T(1)-T0)/Nsim);
gvec(:,1) = exp(  -g_r*sum(log(1-psi_r/c_r))*(T0-t)/Nsim+...
                  -g_t*sum(log(1+(g_l/c_t)*log(1-psi_l/c_l)))*(T0-t)/Nsim+...
                   xi_g  )*P;               

xir_h = xi_r(T(1)-t,r,-1,0,theta_r);
xil_h = xi_l(T(1)-t,lambda,-1,0,theta_l);
xi_h = reshape((xir_h+xil_h')',(N+1)*(L+1),1);
psi_r = psir(T(1),t+U*(T(1)-t),-1,-1,0,0,theta_r,theta_l,rho);
psi_l = psil(T(1),t+U*(T(1)-t),-1,0,theta_l);
hvec(:,1) = exp(  -g_r*sum(log(1-psi_r/c_r))*(T(1)-t)/Nsim+...
                  -g_t*sum(log(1+(g_l/c_t)*log(1-psi_l/c_l)))*(T(1)-t)/Nsim+...
                   xi_h  );
for ell=2:M
    alpha3 = -rr(theta_r,T(ell),T(ell-1));
    xir_g = xi_r(T(ell-1)-t,r,-1,alpha3,theta_r);
    xil_g = xi_l(T(ell-1)-t,lambda,-1,0,theta_l);
    xi_g = reshape((xir_g+xil_g')',(N+1)*(L+1),1);
    psi_r = psir(T(ell-1),t+U*(T(ell-1)-t),-1,-1,alpha3,0,theta_r,theta_l,rho);
    psi_l = psil(T(ell-1),t+U*(T(ell-1)-t),-1,0,theta_l);
    P = exp(  -g_r*sum(log(1+rr(theta_r,T(ell),T(ell-1)+U*(T(ell)-T(ell-1)))/c_r))*...
               (T(ell)-T(ell-1))/Nsim  );
    gvec(:,ell) = exp(  -g_r*sum(log(1-psi_r/c_r))*(T(ell-1)-t)/Nsim+...
                        -g_t*sum(log(1+(g_l/c_t)*log(1-psi_l/c_l)))*(T(ell-1)-t)/Nsim+...
                         xi_g  )*P;

    xir_h = xi_r(T(ell)-t,r,-1,0,theta_r);
    xil_h = xi_l(T(ell)-t,lambda,-1,0,theta_l);
    xi_h = reshape((xir_h+xil_h')',(N+1)*(L+1),1);
    psi_r = psir(T(ell),t+U*(T(ell)-t),-1,-1,0,0,theta_r,theta_l,rho);
    psi_l = psil(T(ell),t+U*(T(ell)-t),-1,0,theta_l);
    hvec(:,ell) = exp(  -g_r*sum(log(1-psi_r/c_r))*(T(ell)-t)/Nsim+...
                        -g_t*sum(log(1+(g_l/c_t)*log(1-psi_l/c_l)))*(T(ell)-t)/Nsim+...
                         xi_h  );
end
TT = T-[T0,T(1:end-1)];
PI = kappa*hvec*TT'-delta*(gvec-hvec)*Id_M;
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
function [llambdau] = llambda(theta_l,t,u)
llambdau = (1-exp(-theta_l*(t-u)))/theta_l;
end