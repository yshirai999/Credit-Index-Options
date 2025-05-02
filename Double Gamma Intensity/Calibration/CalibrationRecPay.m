%% Calibration 01/02/2020
clear
clc

%% Data and Parameters
Treasuryxlsm = fullfile(parentFolder, 'Data', 'TreasuryYield2020.xlsm'); % Adjust 'Data' and file name as needed
CDXOxlsm = fullfile(parentFolder, 'Data', 'IGOptionsData20200102.xlsm'); % Adjust 'Data' and file name as needed
Data_r = xlsread(Treasuryxlsm);
Data_l = xlsread(CDXOxlsm);

% Short rate model guess paramaters
paramin_r = [0.1,15,0.1];
paramax_r = [10,1000,10];
r0_guess = 0;
theta_r_guess = 0.1;
c_r_guess = 400;
g_r_guess = 4;

%% Default intensity model guess parameters
rhomin = 0.1;
paramin_l = [0.01,0,0,0,0];
paramax_l = [30,10000,10000,10000,10000];

theta_l_guess = zeros(1,6);
rho_guess = zeros(1,6);
c_l_guess = zeros(1,6);
g_l_guess = zeros(1,6);
c_t_guess = zeros(1,6);
g_t_guess = zeros(1,6);

% m=1
% ITM
% lambda0 = 1.811898e-02; theta_l = 6.752119e-01; rho = 1.648209e-01; c_l = 4.563077e+02; g_l = 5.364237e+00; c_t = 7.400000e+02; g_t = 5.294346e+00;
% lambda0 = 4.127932e-02; theta_l = 1.538445e+00; rho = 1.611700e-01; c_l = 2.648873e+01; g_l = 5.897348e+00; c_t = 6.939984e+02; g_t = 3.557822e+00;
% lambda0 = 0;            theta_l = 2.959284e+00; rho = 1.078006e-01; c_l = 1.981975e+01; g_l = 7.935774e+00; c_t = 9.304660e+01; g_t = 5.047903e+00;
% theta_l_guess(1) = theta_l; rho_guess(1) = rho-rhomin; c_l_guess(1) = c_l; g_l_guess(1) = g_l; c_t_guess(1) = c_t; g_t_guess(1) = g_t;    
% OTM
lambda0 = 1.801805e-02; theta_l = 7.869028e-01; rho = 1.562220e-01; c_l = 2.032924e+01; g_l = 4.122366e+00; c_t = 6.040000e+02; g_t = 3.319264e+00;
%lambda0 = 1.836509e-02; theta_l = 8.074465e-01; rho = 1.616815e-01; c_l = 1.967187e+01; g_l = 3.977824e+00; c_t = 5.940047e+02; g_t = 3.319211e+00;
theta_l_guess(1) = theta_l; rho_guess(1) = rho-rhomin; c_l_guess(1) = c_l; g_l_guess(1) = g_l; c_t_guess(1) = c_t; g_t_guess(1) = g_t;

% m=2
% ITM
% lambda0 = 4.510659e-03; theta_l = 2.397952e+00; rho = 1.625756e-01; c_l = 3.593098e+00; g_l = 6.858153e+00; c_t = 3.881557e+02; g_t = 3.435002e+00;
% lambda0 = 1.038361e-02; theta_l = 2.969723e+00; rho = 1.404119e-01; c_l = 2.836321e+00; g_l = 7.076869e+00; c_t = 4.279130e+02; g_t = 3.434851e+00;
% theta_l_guess(2) = theta_l; rho_guess(2) = rho-rhomin; c_l_guess(2) = c_l; g_l_guess(2) = g_l; c_t_guess(2) = c_t; g_t_guess(2) = g_t;
% OTM
lambda0 = 0; theta_l = 3.353399e+00; rho = 1.548835e-01; c_l = 4.317856e+00; g_l = 6.061738e+00; c_t = 1.900001e+02; g_t = 3.529832e+00;
% lambda0 = 0; theta_l = 3.398403e+00; rho = 1.223130e-01; c_l = 5.495528e+00; g_l = 6.848726e+00; c_t = 1.653808e+02; g_t = 3.524762e+00;
theta_l_guess(2) = theta_l; rho_guess(2) = rho-rhomin; c_l_guess(2) = c_l; g_l_guess(2) = g_l; c_t_guess(2) = c_t; g_t_guess(2) = g_t;
    
% m=3
%ITM
% lambda0 = 1.912448e-02; theta_l = 2.320101e+00; rho = 1.007303e-01; c_l = 2.913890e+00; g_l = 3.821044e+00; c_t = 3.261490e+02; g_t = 3.382472e+00;
% theta_l_guess(3) = theta_l; rho_guess(3) = rho-rhomin; c_l_guess(3) = c_l; g_l_guess(3) = g_l; c_t_guess(3) = c_t; g_t_guess(3) = g_t;
%OTM
lambda0 = 1.946458e-02; theta_l = 3.229390e+00; rho = 1.118948e-01; c_l = 4.935373e+00; g_l = 2.938284e+00; c_t = 1.095427e+02; g_t = 3.612730e+00;
lambda0 = 1.886404e-02; theta_l = 2.678939e+00; rho = 1.115497e-01; c_l = 6.131311e+00; g_l = 2.698303e+00; c_t = 1.012590e+02; g_t = 3.612362e+00;
% lambda0 = 1.894103e-02; theta_l = 2.770487e+00; rho = 1.104011e-01; c_l = 5.928530e+00; g_l = 3.062303e+00; c_t = 1.139999e+02; g_t = 3.612889e+00;
theta_l_guess(3) = theta_l; rho_guess(3) = rho-rhomin; c_l_guess(3) = c_l; g_l_guess(3) = g_l; c_t_guess(3) = c_t; g_t_guess(3) = g_t;

% m=4
%ITM
% lambda0 = 5.137530e-03; theta_l = 1.233114e+00; rho = 1.000000e-01; c_l = 4.415111e+00; g_l = 3.156418e+00; c_t = 2.574902e+02; g_t = 3.210091e+00;
% lambda0 = 5.068955e-03; theta_l = 1.011975e+00; rho = 1.000000e-01; c_l = 5.627951e+00; g_l = 3.780537e+00; c_t = 2.990158e+02; g_t = 3.210091e+00;
% theta_l_guess(4) = theta_l; rho_guess(4) = rho-rhomin; c_l_guess(4) = c_l; g_l_guess(4) = g_l; c_t_guess(4) = c_t; g_t_guess(4) = g_t;
%OTM
lambda0 = 0; theta_l = 1.793474e-02; rho = 2.651822e-01; c_l = 1.253142e+01; g_l = 3.787548e+00; c_t = 1.286425e+03; g_t = 2.592527e+00;
lambda0 = 0; theta_l = 2.645226e-02; rho = 1.280662e-01; c_l = 1.877568e+01; g_l = 5.183664e+00; c_t = 3.125091e+02; g_t = 2.590363e+00;
% lambda0 = 0; theta_l = 1.000000e-02; rho = 1.294779e-01; c_l = 2.008911e+01; g_l = 5.252971e+00; c_t = 3.114981e+02; g_t = 2.590338e+00;
theta_l_guess(4) = theta_l; rho_guess(4) = rho-rhomin; c_l_guess(4) = c_l; g_l_guess(4) = g_l; c_t_guess(4) = c_t; g_t_guess(4) = g_t;

% m=5
%ITM
% lambda0 = 0; theta_l = 1.000000e-02; rho = 1.053834e-01; c_l = 1.020443e+01; g_l = 3.578843e+00; c_t = 7.770082e+02; g_t = 4.708954e+00;
% lambda0 = 0; theta_l = 4.058194e-02; rho = 1.000211e-01; c_l = 1.053939e+01; g_l = 3.922556e+00; c_t = 7.719489e+02; g_t = 4.836615e+00;
% theta_l_guess(5) = theta_l; rho_guess(5) = rho-rhomin; c_l_guess(5) = c_l; g_l_guess(5) = g_l; c_t_guess(5) = c_t; g_t_guess(5) = g_t;    
%OTM
lambda0 = 0; theta_l = 1.000000e-02; rho = 1.000400e-01; c_l = 1.009812e+01; g_l = 4.420553e+00; c_t = 8.181465e+02; g_t = 4.9855970e+00;
%lambda0 = 0; theta_l = 1.000002e-02; rho = 1.000395e-01; c_l = 1.012016e+01; g_l = 4.642870e+00; c_t = 8.314620e+02; g_t = 4.818457e+00;
theta_l_guess(5) = theta_l; rho_guess(5) = rho-rhomin; c_l_guess(5) = c_l; g_l_guess(5) = g_l; c_t_guess(5) = c_t; g_t_guess(5) = g_t;    

% m=6
%ITM
% lambda0 = 0; theta_l = 3.263195e-01; rho = 1.000001e-01; c_l = 5.956448e+01; g_l = 1.324352e+01; c_t = 2.618867e+02; g_t = 4.308340e+00;
% lambda0 = 0; theta_l = 6.130971e-02; rho = 1.000001e-01; c_l = 4.566360e+01; g_l = 7.477003e+00; c_t = 5.977676e+02; g_t = 8.106857e+00;
% theta_l_guess(6) = theta_l; rho_guess(6) = rho-rhomin; c_l_guess(6) = c_l; g_l_guess(6) = g_l; c_t_guess(6) = c_t; g_t_guess(6) = g_t;
%OTM
% lambda0 = 0; theta_l = 6.130971e-02; rho = 1.000001e-01; c_l = 4.566360e+01; g_l = 7.477003e+00; c_t = 5.977676e+02; g_t = 8.106857e+00;
lambda0 = 0; theta_l = 1.000000e-02; rho = 1.000003e-01; c_l = 8.228926e+01; g_l = 1.024105e+00; c_t = 4.583975e+01; g_t = 8.458455e+00;
% lambda0 = 0; theta_l = 1.000000e-02; rho = 1.000002e-01; c_l = 1.031037e+02; g_l = 1.355065e+00; c_t = 4.951181e+01; g_t = 8.184206e+00;
theta_l_guess(6) = theta_l; rho_guess(6) = rho-rhomin; c_l_guess(6) = c_l; g_l_guess(6) = g_l; c_t_guess(6) = c_t; g_t_guess(6) = g_t;
 
%% Underlying CDX terms
CDXMat = 5; % Underlying CDX maturity (years)
dT = 0.5; % coupons frequency
delta = 0.6; % loss rate
Protection = 0; %if 1, FEP; if 0, no FEP;

%% Short Rate Model Parameter Calibration

% Treasury Parameters
T = Data_r(1,2:end);
Y = Data_r(2,2:end);
t = Data_r(2:end,1);

% Montecarlo integration
Nsim = 100; rng('default'); rng(1); U = rand(1,Nsim);

%  Calibration
w = [1,1,1,1,1,0.5,0.25,0.1,0.05,0.01,0.01]; %weights
ellmax = 5;%max years considered
fun_r = @(x) TreasuryCalib(abs(x(1)),abs(x(2)),abs(x(3)),abs(x(4)),T,Y,U,w,paramin_r,paramax_r,ellmax,0);
x0_r = [r0_guess,theta_r_guess,c_r_guess,g_r_guess];
options = optimset('MaxFunEval',5000,'MaxIter',5000,'Display','off');
[x1_r,fminres_r] = fminunc(fun_r,x0_r,options);

% Short Rate Model Result of Calibration
r0 = abs(x1_r(1));
theta_r = abs(x1_r(2));
c_r = abs(x1_r(3));
g_r = abs(x1_r(4));

range = paramax_r-paramin_r;
if paramin_r(1)>theta_r||theta_r>paramax_r(1)
    n = floor((theta_r-paramin_r(1))/range(1));
    if mod(n,2)==1
        theta_r=theta_r-n*range(1);
    else
        theta_r=paramax_r(1)+n*range(1)-(theta_r-paramin_r(1));
    end
end

if paramin_r(2)>c_r||c_r>paramax_r(2)
    n = floor((c_r-paramin_r(2))/range(2));
    if mod(n,2)==1
        c_r=c_r-n*range(2);
    else
        c_r=paramax_r(2)+n*range(2)-(c_r-paramin_r(2));
    end
end

if paramin_r(3)>g_r||g_r>paramax_r(3)
    n = floor((g_r-paramin_r(3))/range(3));
    if mod(n,2)==1
        g_r=g_r-n*range(3);
    else
        g_r=paramax_r(3)+n*range(3)-(g_r-paramin_r(3));
    end
end

% Short Rate Model Results Visualization
fprintf('\nParameters calibrated\n')
fprintf('r0 = %d; theta_r = %d; c_r = %d; g_r = %d; res = %d\n',r0,theta_r,c_r,g_r,fminres_r)

fprintf('Treasury Market calibration result\n')
TreasuryCalib(r0,theta_r,c_r,g_r,T,Y,U,w,paramin_r,paramax_r,ellmax,1);

%% Default Intensity Model Parameter Calibration

% CDX Parameters
Today = Data_l(1,1);
[Maturity,ind] = unique(Data_l(:,4));
T0 = (Maturity-Today)/360;
M = length(T0);
ind = [ind;length(Data_l(:,3))+1];
CDX_0 = Data_l(1,11);

% Short rate simulation
Nsim = 100;
dt = T0/Nsim;
Gr = ones(Nsim,Nsim,M);
for m = 1:M
    rng('default'); rng(1);
    Gr(:,:,m) = gamrnd(g_r*dt(m),1/c_r,Nsim,Nsim);
end
seed = rng;

%  Calibration
theta_l = zeros(1,M);
rho = zeros(1,M);
c_l = zeros(1,M);
g_l = zeros(1,M);
c_t = zeros(1,M);
g_t = zeros(1,M);
lambda0 = zeros(1,M);
clear 'global'
global ModelPrice_CalibrationRecPay
global MidPrice_CalibrationRecPay
global kappa_CalibrationRecPay
global kappa
global MidRec
global MidPay

Vol = zeros(1,M);
Cub = zeros(1,M);
Qua = zeros(1,M);
Vol_Model = zeros(1,M);
Cub_Model = zeros(1,M);
Qua_Model = zeros(1,M);
Vol_Model_Check = zeros(1,M);
Cub_Model_Check = zeros(1,M);
Qua_Model_Check = zeros(1,M);

Var = zeros(1,M);
Skew = zeros(1,M);
Kurt = zeros(1,M);
Var_Model = zeros(1,M);
Skew_Model = zeros(1,M);
Kurt_Model = zeros(1,M);
Var_Model_Check = zeros(1,M);
Skew_Model_Check = zeros(1,M);
Kurt_Model_Check = zeros(1,M);

%
for m=1:6
    Gr_m = Gr(:,:,m);
    fun_l = @(x) CDXOMCCalibRecPay(r0,delta,U,Gr_m,seed,...
        theta_r,abs(x(1)),abs(x(2))+rhomin,c_r,abs(x(3)),abs(x(5)),g_r,abs(x(4)),...
        abs(x(6)),Data_l,T0,ind,dT,CDXMat,m,CDX_0,paramin_l,paramax_l,2,...
        Protection);
    x0_l = [theta_l_guess(m),rho_guess(m),c_l_guess(m),g_l_guess(m),...
        c_t_guess(m),g_t_guess(m)];
    options = optimset('MaxFunEval',5000,'MaxIter',5000,'Display','off');
    %[x1_l,fminres_l] = fminsearch(fun_l,x0_l,options);
    x1_l = x0_l;
    
    % Default Intensity Model Results of Calibration
    
    theta_l(m) = abs(x1_l(1));
    rho(m) = abs(x1_l(2))+rhomin;
    c_l(m) = abs(x1_l(3));
    g_l(m) = abs(x1_l(4));
    c_t(m) = abs(x1_l(5));
    g_t(m) = abs(x1_l(6));
    range = paramax_l-paramin_l;
    
    if paramin_l(1)>theta_l(m)||theta_l(m)>paramax_l(1)
        n = floor((theta_l(m)-paramin_l(1))/range(1));
        if mod(n,2)==1
            theta_l(m)=theta_l(m)-n*range(1);
        else
            theta_l(m)=paramax_l(1)+n*range(1)-(theta_l(m)-paramin_l(1));
        end
    end
    if paramin_l(2)>c_l(m)||c_l(m)>paramax_l(2)
        n = floor((c_l(m)-paramin_l(2))/range(2));
        if mod(n,2)==1
            c_l(m)=c_l(m)-n*range(2);
        else
            c_l(m)=paramax_l(2)+n*range(2)-(c_l(m)-paramin_l(2));
        end
    end
    if paramin_l(3)>g_l(m)||g_l(m)>paramax_l(3)
        n = floor((g_l(m)-paramin_l(3))/range(3));
        if mod(n,2)==1
            g_l(m)=g_l(m)-n*range(3);
        else
            g_l(m)=paramax_l(3)+n*range(3)-(g_l(m)-paramin_l(3));
        end
    end
    if paramin_l(4)>c_t(m)||c_t(m)>paramax_l(4)
        n = floor((c_t(m)-paramin_l(4))/range(4));
        if mod(n,2)==1
            c_t(m)=c_t(m)-n*range(4);
        else
            c_t(m)=paramax_l(4)+n*range(4)-(c_t(m)-paramin_l(4));
        end
    end
    if paramin_l(5)>g_t(m)||g_t(m)>paramax_l(5)
        n = floor((g_t(m)-paramin_l(5))/range(5));
        if mod(n,2)==1
            g_t(m)=g_t(m)-n*range(5);
        else
            g_t(m)=paramax_l(5)+n*range(5)-(g_t(m)-paramin_l(5));
        end
    end
    
    % Results Visualization
    options = optimset('Display','off');
    lambda0(m) = fsolve(@(x)CDXPIDE(0,0,r0,abs(x),[dT:dT:CDXMat],...
         CDX_0*1e-4,delta,U,theta_r,theta_l(m),rho(m),c_r,c_l(m),c_t(m),...
         g_r,g_l(m),g_t(m),0,0,0),0,options);
    fprintf('lambda0 = %d; theta_l = %d; rho = %d; c_l = %d; g_l = %d; c_t = %d; g_t = %d;\n',lambda0(m),theta_l(m),rho(m),c_l(m),g_l(m),c_t(m),g_t(m))
    fprintf('CDX Market calibration result\n')
    CDXOMCCalibRecPay(r0,delta,U,Gr_m,seed,...
        theta_r,theta_l(m),rho(m),c_r,c_l(m),c_t(m),g_r,g_l(m),g_t(m),Data_l,T0,ind,dT,CDXMat,m,CDX_0,paramin_l,paramax_l,1,...
        Protection); 
    
    % Market Implied Annuity Price and Forward Spread
    payer_m_mkt = MidPay*1e+4;
    receiver_m_mkt = MidRec*1e+4;
    strike_m = kappa*1e+4;
    c_hat = strike_m(end);
    c_f_mkt = fsolve(@(x)interp1(strike_m,payer_m_mkt-receiver_m_mkt,x),CDX_0,options);
    A_0_mkt = interp1(strike_m,payer_m_mkt-receiver_m_mkt,0,'linear','extrap')/c_f_mkt;
    cmin = min(strike_m);
    cmax = max(strike_m);
    
    [~,kappa_f_ind]=min((strike_m-c_f_mkt).^2); %Assume chat = c_f
    if strike_m(kappa_f_ind)>c_f_mkt
        strike_m_receiver = [strike_m(1:kappa_f_ind-1);c_f_mkt];
        strike_m_payer = [c_f_mkt;strike_m(kappa_f_ind:end)];
    else
        strike_m_receiver = [strike_m(1:kappa_f_ind);c_f_mkt];
        strike_m_payer = [c_f_mkt;strike_m(kappa_f_ind+1:end)];
    end
    
    receiver_m_mkt = interp1(strike_m,receiver_m_mkt,strike_m_receiver,'linear','extrap');
    payer_m_mkt = interp1(strike_m,payer_m_mkt,strike_m_payer,'linear','extrap');
    
    Vol(m) = trapz(strike_m_receiver,2*receiver_m_mkt)...
              +trapz(strike_m_payer,2*payer_m_mkt);
    Cub(m) = trapz(strike_m_receiver,6*(strike_m_receiver-c_f_mkt).*receiver_m_mkt)...
              +trapz(strike_m_payer,6*(strike_m_payer-c_f_mkt).*payer_m_mkt);
    Qua(m) = trapz(strike_m_receiver,(12*(strike_m_receiver-c_f_mkt).^2).*receiver_m_mkt)...
              +trapz(strike_m_payer,(12*(strike_m_payer-c_f_mkt).^2).*payer_m_mkt);
        
	Var(m) = Vol(m)/A_0_mkt;
    Skew(m) = ( Cub(m)/A_0_mkt ) / ( Var(m)^(3/2) );
    Kurt(m) = ( Qua(m)/A_0_mkt ) / ( Var(m)^(2) );
    
    %Model Implied Option Price Surface
    rng('default'); rng(seed);
    Gt_m = gamrnd(g_t(m)*dt(m),1/c_t(m),Nsim,Nsim);
    Gl_m = gamrnd(g_l(m)*Gt_m,1/c_l(m),Nsim,Nsim);
    er = exp( -theta_r * ( T0(m)-cumsum(ones(Nsim,1))*dt(m) ) );
    el = exp( -theta_l(m) * ( T0(m)-cumsum(ones(Nsim,1))*dt(m) ) );
    eY_r = ( 1 - exp( -theta_r * ( T0(m)-cumsum(ones(Nsim,1))*dt(m) ) ) ) / theta_r;
    rT0 = r0 * exp ( -theta_r*T0(m) ) + Gr_m*er; %simulated short rate @ T0(m)
    lT0 = lambda0(m) * exp( -theta_l(m)*T0(m) ) + (rho(m)*Gr_m+Gl_m)*el; %simulated default intensity @ T0(m)
    Y_rT0 = r0 * ( 1-exp(-theta_r*T0(m)) ) / theta_r + Gr_m*eY_r; %simulated integrated short rate @ T0(m)
    
    strike_m_mdl = kappa;
    payer_m_mdl = zeros(length(strike_m_mdl),1);
    receiver_m_mdl = zeros(length(strike_m_mdl),1);
    
    for kk = 1:length(strike_m_mdl)
        for s=1:Nsim
            payer_m_mdl(kk) = payer_m_mdl(kk) +  exp(-Y_rT0(s))...
                *max(-CDXPIDE(T0(m),T0(m),rT0(s),lT0(s),[T0(m)+dT:dT:CDXMat],strike_m_mdl(kk),delta,U,...
                theta_r,theta_l(m),rho(m),c_r,c_l(m),c_t(m),g_r,g_l(m),g_t(m),0,0,0),0)/Nsim;
            receiver_m_mdl(kk) = receiver_m_mdl(kk) +  exp(-Y_rT0(s))...
                *max(CDXPIDE(T0(m),T0(m),rT0(s),lT0(s),[T0(m)+dT:dT:CDXMat],strike_m_mdl(kk),delta,U,...
                theta_r,theta_l(m),rho(m),c_r,c_l(m),c_t(m),g_r,g_l(m),g_t(m),0,0,0),0)/Nsim;
        end
        fprintf('strike = %d, payer = %d, receiver = %d\n',strike_m_mdl(kk),payer_m_mdl(kk),receiver_m_mdl(kk))
    end
    
    % Model Implied Annuity Price
    payer_m_mdl = payer_m_mdl*1e+4;
    receiver_m_mdl = receiver_m_mdl*1e+4;
    strike_m_mdl = strike_m_mdl*1e+4;
    c_f_mdl = fsolve(@(x)interp1(strike_m_mdl,payer_m_mdl-receiver_m_mdl,x,'linear','extrap'),CDX_0,options);
    A_0_mdl = interp1(strike_m_mdl,payer_m_mdl-receiver_m_mdl,0,'linear','extrap')/c_f_mdl;
    
    % Model Implied Spread Moments
    
    receiver_m_mdl = interp1(strike_m_mdl,receiver_m_mdl,strike_m_receiver,'linear','extrap');
    payer_m_mdl = interp1(strike_m_mdl,payer_m_mdl,strike_m_payer,'linear','extrap');
    
    Vol_Model(m) = trapz(strike_m_receiver,2*receiver_m_mdl)...
              +trapz(strike_m_payer,2*payer_m_mdl);
    Cub_Model(m) = trapz(strike_m_receiver,6*(strike_m_receiver-c_f_mdl).*receiver_m_mdl)...
              +trapz(strike_m_payer,6*(strike_m_payer-c_f_mdl).*payer_m_mdl);
    Qua_Model(m) = trapz(strike_m_receiver,(12*(strike_m_receiver-c_f_mdl).^2).*receiver_m_mdl)...
              +trapz(strike_m_payer,(12*(strike_m_payer-c_f_mdl).^2).*payer_m_mdl);
    
    Var_Model(m) = Vol_Model(m)/A_0_mdl;
    Skew_Model(m) = ( Cub_Model(m)/A_0_mdl ) / ( Var_Model(m)^(3/2) );
    Kurt_Model(m) = ( Qua_Model(m)/A_0_mdl ) / ( Var_Model(m)^(2) );
    
    % Check
    A_0_mdl_Check = Annuity(0,T0(m),r0,lambda0(m),[T0(m)+dT:dT:CDXMat],U,...
            theta_r,theta_l(m),rho(m),c_r,c_l(m),g_r,g_l(m),0,0);
    fp_m_0 = -CDXPIDE(0,T0(m),r0,lambda0(m),[T0(m)+dT:dT:CDXMat],0,delta,U,...
                theta_r,theta_l(m),rho(m),c_r,c_l(m),c_t(m),g_r,g_l(m),g_t(m),0,0,Protection);
    c_f_mdl_Check = (fp_m_0/A_0_mdl_Check)*1e+4;
    P=0;
    for s=1:Nsim
        c_T0 = CDXSpread(T0(m),T0(m),rT0(s),lT0(s),[T0(m)+dT:dT:CDXMat],delta,U,...
            theta_r,theta_l(m),rho(m),c_r,c_l(m),c_t(m),g_r,g_l(m),g_t(m),0,0,0)*1e+4; %No FEP for spot index (at time T0(m))
        A = Annuity(T0(m),T0(m),rT0(s),lT0(s),[T0(m)+dT:dT:CDXMat],U,...
            theta_r,theta_l(m),rho(m),c_r,c_l(m),g_r,g_l(m),0,0);
        Vol_Model_Check(m) = Vol_Model_Check(m) + exp(-Y_rT0(s))*(A/A_0_mdl_Check)*((c_T0-c_f_mdl_Check)^2*(c_T0<max(strike_m_mdl))*(c_T0>min(strike_m_mdl)))/Nsim;
        Cub_Model_Check(m) = Cub_Model_Check(m) + exp(-Y_rT0(s))*(A/A_0_mdl_Check)*((c_T0-c_f_mdl_Check)^3*(c_T0<max(strike_m_mdl))*(c_T0>min(strike_m_mdl)))/Nsim;
        Qua_Model_Check(m) = Qua_Model_Check(m) + exp(-Y_rT0(s))*(A/A_0_mdl_Check)*((c_T0-c_f_mdl_Check)^4*(c_T0<max(strike_m_mdl))*(c_T0>min(strike_m_mdl)))/Nsim;
        P = P+(c_T0>c_f_mdl_Check)/Nsim;
    end
    
    Vol_Model_Check(m) = A_0_mdl_Check * ( Vol_Model_Check(m));
    Cub_Model_Check(m) = A_0_mdl_Check * ( Cub_Model_Check(m));
    Qua_Model_Check(m) = A_0_mdl_Check * ( Qua_Model_Check(m));
    
    Var_Model_Check(m) = Vol_Model_Check(m)/A_0_mdl_Check;
    Skew_Model_Check(m) = ( Cub_Model_Check(m)/A_0_mdl_Check ) / ( Var_Model_Check(m)^(3/2) );
    Kurt_Model_Check(m) = Qua_Model_Check(m)/A_0_mdl_Check / ( Var_Model_Check(m)^(2) );

end

%% Scatter Plot

[K,~] = size((ModelPrice_CalibrationRecPay));
kappa = zeros(K,1);
kappa(1:K)=kappa_CalibrationRecPay(1:K);

% Specify a folder for saving plots
% Get the current folder
currentFolder = pwd;
% Get the parent folder
[parentFolder, ~, ~] = fileparts(currentFolder);
% Specify the target folder in the parent directory
fpath = fullfile(parentFolder, 'Plots');
% Create the folder if it doesn't exist  
if ~exist(fpath, 'dir')
    mkdir(fpath);
end 

figure
hold on
box on
grid on
Colors = {'k','b','r','g',[.4 .1 .6],[.5 .6 .7],[.8 .2 .6]}; % Cell array of colors.
for k=1:K
    model = (ModelPrice_CalibrationRecPay(k,:)*1e-4)';
    market = (MidPrice_CalibrationRecPay(k,:)*1e-4)';
    kappa_k = kappa(k)*ones(1,length(T0));
    s = scatter3(kappa_k(market>0),T0(market>0),model(market>0),'MarkerEdgeColor',Colors{k});
    s.Annotation.LegendInformation.IconDisplayStyle = 'off';
    scatter3(kappa_k(market>0),T0(market>0),market(market>0),'*','MarkerEdgeColor',Colors{k});
    %str_k{k} = strcat('Strike = ',num2str(kappa(k)));
end
%legend (str_k,'interpreter','latex','location','southeast');    
view(6,12)
xlabel('Strike','interpreter','latex')
ylabel('Term','interpreter','latex')
hold off
str=strcat('CalibrationResult');
fname=str;
saveas(gcf, fullfile(fpath, fname), 'epsc');
close all

%% Variance, Skewness Kurtosis
fprintf('\n Var, Skewness, Kurtosis comparison\n')
for m = 1:M
    fprintf('Variance Mkt = %d, Variance Model = %d, Variance Model MC = %d\n', Var(m), Var_Model(m), Var_Model_Check(m));
end
for m = 1:M
    fprintf('Skeness Mkt = %d, Skeness Model = %d, Skeness Model MC = %d\n', Skew(m), Skew_Model(m), Skew_Model_Check(m));
end
for m = 1:M
    fprintf('Kurtosis Mkt = %d, Kurtosis Model = %d, Kurtosis Model MC = %d\n', Kurt(m), Kurt_Model(m), Kurt_Model_Check(m));
end
