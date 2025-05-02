function [R] = Rfun(r,lambda,uu,Y_r,Y_t,w,Nsim,rho,g_r,g_t,c_r,c_t,N,L)
%Compute integrals in RHS of PIDE. Input u is the CDXO price surface at t.

% Jumps in default intensity
R_l = zeros(L-1,N-1);
for k=2:L
    LL = lambda(k)+Y_t';
    u_l = interp1(lambda(k:L+1)',uu(k:L+1,2:N),LL,'linear','extrap');
    R_l(k-1,:) = sum( ((u_l-uu(k,2:N)).*w)./Y_t' );
end
R_l = (g_t/c_t)*R_l/Nsim;

% Jumps in interest rates
R_r = zeros(L-1,N-1);
for d=N-2:-1:-N+2 %we construct each diagonal of matrix R_r
    imin = max(d+1,1); %minimum value of r along diagonal d
    imax = min(N+1+d,N+1); %max value of r along diagonal d
    ud = diag(uu,d);
    RR = r(imin+1:imax-1)+Y_r'; %i-th columns are jumps for starting point r_i
    u_r = max(interp1(sqrt(1+rho^2)*r(imin:imax)',ud,sqrt(1+rho^2)*RR,'linear','extrap'),0); % matrix of size as RR
    dd = sum((u_r-ud(2:end-1)')./Y_r');
    R_r = R_r+diag(dd,d);
end
R_r = (g_r/c_r)*R_r/Nsim;

% Integral
R = reshape((R_r)',(N-1)*(L-1),1);
R = R + reshape((R_l)',(N-1)*(L-1),1);
end

