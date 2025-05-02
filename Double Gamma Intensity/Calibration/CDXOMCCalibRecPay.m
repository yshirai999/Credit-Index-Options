function [err] = CDXOMCCalibRecPay(r,delta,U,Gr_m,seed,...
    theta_r,theta_l,rho,c_r,c_l,c_t,g_r,g_l,g_t,Data,T0,ind,dT,CDXMat,m,CDX_0,...
    paramin,paramax,Display,Protection)
%r(T_0) and lambda(T_0) simulated via MonteCarlo, and CDX payoff computed
%accordingly.

global ModelPrice_CalibrationRecPay
global MidPrice_CalibrationRecPay
global AskPrice_CalibrationRecPay
global BidPrice_CalibrationRecPay
global kappa_CalibrationRecPay

global MidRec
global MidPay
global kappa

range = paramax-paramin;

if paramin(1)>theta_l||theta_l>paramax(1)
    n = floor((theta_l-paramin(1))/range(1));
    if mod(n,2)==1
        theta_l=theta_l-n*range(1);
    else
        theta_l=paramax(1)+n*range(1)-(theta_l-paramin(1));
    end
end

if paramin(2)>c_l||c_l>paramax(2)
    n = floor((c_l-paramin(2))/range(2));
    if mod(n,2)==1
        c_l=c_l-n*range(2);
    else
        c_l=paramax(2)+n*range(2)-(c_l-paramin(2));
    end
end

if paramin(3)>g_l||g_l>paramax(3)
    n = floor((g_l-paramin(3))/range(3));
    if mod(n,2)==1
        g_l=g_l-n*range(3);
    else
        g_l=paramax(3)+n*range(3)-(g_l-paramin(3));
    end
end

if paramin(4)>c_t||c_t>paramax(4)
    n = floor((c_t-paramin(4))/range(4));
    if mod(n,2)==1
        c_t=c_t-n*range(4);
    else
        c_t=paramax(4)+n*range(4)-(c_t-paramin(4));
    end
end

if paramin(5)>g_t||g_t>paramax(5)
    n = floor((g_t-paramin(5))/range(5));
    if mod(n,2)==1
        g_t=g_t-n*range(5);
    else
        g_t=paramax_l(5)+n*range(5)-(g_t-paramin_l(5));
    end
end

options = optimset('Display','off');
lambda = abs( fsolve(@(x)CDXPIDE(0,0,r,abs(x),[dT:dT:CDXMat],CDX_0*1e-4,delta,U,...
                theta_r,theta_l,rho,c_r,c_l,c_t,g_r,g_l,g_t,0,0,0),0,options) ); %Note: no front end protection for spot CDX index.
%fprintf('lambda0 = %d; theta_l = %d; rho = %d; c_l = %d; g_l = %d;\n',lambda,theta_l,rho,c_l,g_l)

% CDX data
err = 0;

T = T0(m)+dT:dT:CDXMat;
kappa = Data(ind(m):ind(m+1)-1,3)*1e-4;
MidPay = Data(ind(m):ind(m+1)-1,6)*1e-4;
BidPay = Data(ind(m):ind(m+1)-1,7)*1e-4;
AskPay = Data(ind(m):ind(m+1)-1,8)*1e-4;
MidRec = Data(ind(m):ind(m+1)-1,18)*1e-4;
BidRec = Data(ind(m):ind(m+1)-1,19)*1e-4;
AskRec = Data(ind(m):ind(m+1)-1,20)*1e-4;

% MC Simulation
Nsim = length(Gr_m);
dt = T0(m)/Nsim;
rng('default'); rng(seed);
Gt_m = gamrnd(g_t*dt,1/c_t,[Nsim,Nsim]);
Gl_m = gamrnd(g_l*Gt_m,1/c_l,[Nsim,Nsim]);
er = exp( -theta_r * ( T0(m)-cumsum(ones(Nsim,1))*dt ) );
el = exp( -theta_l * ( T0(m)-cumsum(ones(Nsim,1))*dt ) );
eY_r = ( 1 - exp( -theta_r * ( T0(m)-cumsum(ones(Nsim,1))*dt ) ) ) / theta_r;
rT0 = r * exp ( -theta_r*T0(m) ) + Gr_m*er; %simulated short rate @ T0(m)
lT0 = lambda * exp( -theta_l*T0(m) ) + (rho*Gr_m+Gl_m)*el; %simulated default intensity @ T0(m)
Y_rT0 = r * ( 1-exp(-theta_r*T0(m)) ) / theta_r + Gr_m*eY_r; %simulated int short rate @ T0(m)
    
% Compute Error
for kk = 1:length(kappa)
    Pkk = 0;
    if (abs(kappa(kk)-CDX_0*1e-4)/(CDX_0*1e-4)<0.3)&&(abs(kappa(kk)-CDX_0*1e-4)/(CDX_0*1e-4)>0)
        kappa_CalibrationRecPay(kk) = kappa(kk);
        if kappa(kk)<CDX_0*1e-4  % If <, calibrating to ITM with respect to bond (or, in terms of spread, calibrating to OTM)
            for s=1:Nsim
                Pkk = Pkk + exp(-Y_rT0(s))...
                    *max(CDXPIDE(T0(m),T0(m),rT0(s),lT0(s),T,kappa(kk),delta,U,... 
                    theta_r,theta_l,rho,c_r,c_l,c_t,g_r,g_l,g_t,0,0,0),0)/Nsim; % This is the receiver price
            end
            %fprintf('Receiver: Model Price= %d, Market Price = %d, Strike = %4.2f, g_l = %d, c_l = %d\n', Pkk*1e+4, MidRec(kk)*1e+04, kappa(kk)*1e+04, g_l, c_l)
            if Protection == 1
                Pkk = max(Pkk-EFEP(0,T0(m),r,lambda,delta,U,theta_r,theta_l,rho,c_r,c_l,g_r,g_l,0,0),0);
            end 
            err = err + ( Pkk-MidRec(kk) )^2; %Quadratic error for strike kappa(kk) and maturity T0(m)
            if Display == 1
                fprintf('Type = Receiver, Strike = %4.2f, Term = %1.2f:\n    Mid = %4.2f, Bid = %4.2f, Ask = %4.2f, Model Price = %4.2f\n',...
                    kappa(kk)*1e+04,T0(m), MidRec(kk)*1e+04, BidRec(kk)*1e+04,AskRec(kk)*1e+04,Pkk*1e+04);
            end
            ModelPrice_CalibrationRecPay(kk,m) = Pkk*1e+04;
            MidPrice_CalibrationRecPay(kk,m) = MidRec(kk)*1e+04;
            AskPrice_CalibrationRecPay(kk,m) = AskRec(kk)*1e+04;
            BidPrice_CalibrationRecPay(kk,m) = BidRec(kk)*1e+04;
        else  %Payer is OTM, Receiver is ITM
            for s=1:Nsim
                Pkk = Pkk +  exp(-Y_rT0(s))...
                    *max(-CDXPIDE(T0(m),T0(m),rT0(s),lT0(s),T,kappa(kk),delta,U,...
                    theta_r,theta_l,rho,c_r,c_l,c_t,g_r,g_l,g_t,0,0,0),0)/Nsim;
            end
            %fprintf('Payer: Model Price= %d, Market Price = %d, Strike = %4.2f, g_l = %d, c_l = %d\n', Pkk*1e+4, MidPay(kk)*1e+04, kappa(kk)*1e+04, g_l, c_l)
            if Protection == 1
                Pkk = Pkk+EFEP(0,T0(m),r,lambda,delta,U,theta_r,theta_l,rho,c_r,c_l,g_r,g_l,0,0);
            end
            err = err + ( Pkk-MidPay(kk) )^2; %Quadratic error for strike kappa(kk) and maturity T0(m)
            if Display == 1
                fprintf('Type = payer, Strike = %4.2f, Term = %1.2f:\n    Mid = %4.2f, Bid = %4.2f, Ask = %4.2f, Model Price = %4.2f\n',...
                    kappa(kk)*1e+04,T0(m), MidPay(kk)*1e+04, BidPay(kk)*1e+04,AskPay(kk)*1e+04,Pkk*1e+04);
            end
            ModelPrice_CalibrationRecPay(kk,m)=Pkk*1e+04;
            MidPrice_CalibrationRecPay(kk,m)=MidPay(kk)*1e+04;
            AskPrice_CalibrationRecPay(kk,m)=AskPay(kk)*1e+04;
            BidPrice_CalibrationRecPay(kk,m)=BidPay(kk)*1e+04;
        end
    end
end

if Display == 2
    sqrt(err)
end
end

