%% Calibration IG33 to 1-Month options, 01/02/2020-06/05/2020
clear
clc

CalibrationMaturity = 1; % specify maturity of options used to calibrate model

%% Data and Parameters
Treasuryxlsm = fullfile(parentFolder, 'Data', 'TreasuryYield2020.xlsm'); % Adjust 'Data' and file name as needed
CDXOxlsmP = fullfile(parentFolder, 'Data', 'PayerIG332020.xlsm'); % Adjust 'Data' and file name as needed
CDXOxlsmR = fullfile(parentFolder, 'Data', 'ReceiverIG332020.xlsm'); % Adjust 'Data' and file name as needed
Data_r = xlsread(Treasuryxlsm);
Data_P = xlsread(CDXOxlsmP);
Data_R = xlsread(CDXOxlsmR);

% Short rate model paramaters range
paramin_r = [0.1,15,0.1];
paramax_r = [10,1000,10];

% Default intensity model parameters range
rhomin = 0.1;
paramin_l = [0.01,0,0,0,0];
paramax_l = [30,10000,10000,10000,10000];

% Underlying CDX terms
CDXMat = 5; % Underlying CDX maturity (years)
dT = 0.5; % coupons frequency
delta = 0.4; % recovery rate
Protection = 0;

%% Data and initial parameters

% Treasury Parameters
T = Data_r(1,2:end);
t0 = Data_r(2:end,1);
[t0_l,ind_t0_l_P] = unique(Data_P(1:end,1));
[~,ind_t0_l_R] = unique(Data_R(1:end,1));
ind_t0_l_P = [ind_t0_l_P;length(Data_P)+1];
ind_t0_l_R = [ind_t0_l_R;length(Data_R)+1];

% Montecarlo integration
Nsim = 100;
rng('default'); %freeze the seed
rng(1);
U = rand(1,Nsim);

% Variables Definition
Vol_Model = zeros(1,length(t0_l));
Vol = zeros(1,length(t0_l));
Cub_Model = zeros(1,length(t0_l));
Cub = zeros(1,length(t0_l));
Qua_Model = zeros(1,length(t0_l));
Qua = zeros(1,length(t0_l));

c_f_mdl = zeros(1,length(t0_l));
c_f_mkt = zeros(1,length(t0_l));
A_0_mdl = zeros(1,length(t0_l));
A_0_mkt = zeros(1,length(t0_l));

c_hat = zeros(1,length(t0_l));

Var_Model = zeros(1,length(t0_l));
Var = zeros(1,length(t0_l));
Skew_Model = zeros(1,length(t0_l));
Skew = zeros(1,length(t0_l));
Kurt_Model = zeros(1,length(t0_l));
Kurt = zeros(1,length(t0_l));

r0 = zeros(1,length(t0_l));
theta_r = zeros(1,length(t0_l));
c_r = zeros(1,length(t0_l));
g_r = zeros(1,length(t0_l));
lambda0 = zeros(1,length(t0_l));
theta_l = zeros(1,length(t0_l));
rho = zeros(1,length(t0_l));
c_l = zeros(1,length(t0_l));
g_l = zeros(1,length(t0_l));
c_t = zeros(1,length(t0_l));
g_t = zeros(1,length(t0_l));

r0(1) = 0.01468981; theta_r(1) = 0.5500511; c_r(1) = 400.0005; g_r(1) = 3.947560;
lambda0(1) = 0; theta_l(1) = 3.353399e+00; rho(1) = 1.548835e-01;
    c_l(1) = 4.317856e+00; g_l(1) = 6.061738e+00; c_t(1) = 1.900001e+02; g_t(1) = 3.529832e+00;
    
%%  Calibration
t_r = find(t0==t0_l(1));

fprintf('\nParameters calibrated for t0 = %d (%d/%d/%d)\n',t0(t_r),month(t0(t_r)-1),day(t0(t_r)-1),year(t0(t_r)-1)+1900)
fprintf('r0 = %d; theta_r = %d; c_r = %d; g_r = %d;\n',r0(1),theta_r(1),c_r(1),g_r(1))
fprintf('lambda0 = %d; theta_l = %d; rho = %d; c_l = %d; g_l = %d;\n',lambda0(1),theta_l(1),rho(1),c_l(1),g_l(1))

for t=2:length(t0_l)
    
    % Global parameters for no arbitrage implied vol cub quartic contracts
    clear 'global'
    global kappa
    global MidRec
    global MidPay
    
    t_r = find(t0==t0_l(t));
    Y = Data_r(t_r+1,2:end);
    
    w = [1,1,1,1,1,0.5,0.25,0.1,0.05,0.01,0.01]; %weights
    ellmax = 5;%max years considered
    fun_r = @(x) TreasuryCalib(abs(x(1)),abs(x(2)),abs(x(3)),abs(x(4)),T,Y,U,w,paramin_r,paramax_r,ellmax,0);
    x0_r = [r0(t-1),theta_r(t-1),c_r(t-1),g_r(t-1)];
    options = optimset('MaxFunEval',5000,'MaxIter',5000,'Display','off');
    [x1_r,fminres_r] = fminsearch(fun_r,x0_r,options);

    % Short Rate Model Result of Calibration
    r0(t) = abs(x1_r(1));
    theta_r(t) = abs(x1_r(2));
    c_r(t) = abs(x1_r(3));
    g_r(t) = abs(x1_r(4));
    
    range = paramax_r-paramin_r;
    if paramin_r(1)>theta_r(t)||theta_r(t)>paramax_r(1)
        n = floor((theta_r(t)-paramin_r(1))/range(1));
        if mod(n,2)==1
            theta_r(t)=theta_r(t)-n*range(1);
        else
            theta_r(t)=paramax_r(1)+n*range(1)-(theta_r(t)-paramin_r(1));
        end
    end
    
    if paramin_r(2)>c_r(t)||c_r(t)>paramax_r(2)
        n = floor((c_r(t)-paramin_r(2))/range(2));
        if mod(n,2)==1
            c_r(t)=c_r(t)-n*range(2);
        else
            c_r(t)=paramax_r(2)+n*range(2)-(c_r(t)-paramin_r(2));
        end
    end
    
    if paramin_r(3)>g_r(t)||g_r(t)>paramax_r(3)
        n = floor((g_r(t)-paramin_r(3))/range(3));
        if mod(n,2)==1
        g_r(t)=g_r(t)-n*range(3);
            else
        g_r(t)=paramax_r(3)+n*range(3)-(g_r(t)-paramin_r(3));
        end
    end
    
    % Short Rate Model Results Visualization
    fprintf('\nParameters calibrated for t0 = %d (%d/%d/%d)\n',t0(t_r),month(t0(t_r)-1),day(t0(t_r)-1),year(t0(t_r)-1)+1900)
    fprintf('r0 = %d; theta_r = %d; c_r = %d; g_r = %d; res = %d\n',r0(t),theta_r(t),c_r(t),g_r(t),fminres_r)
    fprintf('Treasury Market calibration result\n')
    TreasuryCalib(r0(t),theta_r(t),c_r(t),g_r(t),T,Y,U,w,paramin_r,paramax_r,ellmax,1);
    
    % CDX Parameters
    Today = t0(t_r);
    [AllTerms,indAllTerms] = unique(Data_P(ind_t0_l_P(t):ind_t0_l_P(t+1)-1,6)); %select term (in month)
    if isempty(AllTerms(AllTerms==CalibrationMaturity))
        break
    else
        indMat = indAllTerms(AllTerms==CalibrationMaturity);
    end
    Maturity=Data_P(ind_t0_l_P(t)+indMat-1,5); %Obtain maturity corresponding to term
    T0 = (Maturity-Today)/360;  
    CDX_0 = Data_P(ind_t0_l_P(t),13);
    
    Data_t_P = Data_P(ind_t0_l_P(t):ind_t0_l_P(t+1)-1,:);
    Data_t_R = Data_R(ind_t0_l_R(t):ind_t0_l_R(t+1)-1,:);
    ind_t_P = find(Data_t_P(:,5)==Maturity);
    ind_t_R = find(Data_t_R(:,5)==Maturity);
    Data_t_P = Data_t_P(ind_t_P,:);
    Data_t_R = Data_t_R(ind_t_R,:);    
    Data_t_P = reshape(Data_t_P(~isnan(Data_t_P)),size(Data_t_P)-[0,length(Data_t_P(1,isnan(Data_t_P(1,:))))]);
    Data_t_R = reshape(Data_t_R(~isnan(Data_t_R)),size(Data_t_R)-[0,length(Data_t_R(1,isnan(Data_t_R(1,:))))]);
    
    strikerange = max(Data_t_P(1,3),Data_t_R(1,3)):2.5:min(Data_t_P(end,3),Data_t_R(end,3));
    indstrike_P = find(Data_t_P(:,3)==strikerange(1)):1:find(Data_t_P(:,3)==strikerange(end));
    indstrike_R = find(Data_t_R(:,3)==strikerange(1)):1:find(Data_t_R(:,3)==strikerange(end));
    Data_t_P = Data_t_P(indstrike_P,:);
    Data_t_R = Data_t_R(indstrike_R,:);
    Data_t = [Data_t_P,NaN(size(Data_t_P(:,1))),Data_t_R];
    
    % Short rate simulation
    Nsim = 100;
    dt = T0/Nsim;
    rng('default'); rng(1);
    Gr = gamrnd(g_r(t)*dt,1/c_r(t),[Nsim,Nsim]);
    seed = rng;
    
    % Default Intensity Model Calibration
    ind = [1,length(Data_t_P(:,3))+1];
    format long
    fun_l = @(x) CDXOMCCalibRecPay(r0(t),delta,U,Gr,seed,...
        theta_r(t),abs(x(1)),abs(x(2))+rhomin,c_r(t),abs(x(3)),abs(x(5)),g_r(t),abs(x(4)),...
        abs(x(6)),Data_t,T0,ind,dT,CDXMat,1,CDX_0,paramin_l,paramax_l,0,0);
    x0_l = [theta_l(t-1),rho(t-1),c_l(t-1),g_l(t-1),c_t(t-1),g_t(t-1)];
    options = optimset('MaxFunEval',5000,'MaxIter',100,'Display','off','TolFun',1e-4,'TolX',1e-4);
    [x1_l,fminres_l] = fminsearch(fun_l,x0_l,options);
    
    % Default Intensity Model Calibration Results    
    theta_l(t) = abs(x1_l(1));
    rho(t) = abs(x1_l(2))+rhomin;
    c_l(t) = abs(x1_l(3));
    g_l(t) = abs(x1_l(4));
    c_t(t) = abs(x1_l(5));
    g_t(t) = abs(x1_l(6));
    range = paramax_l-paramin_l;
    
    if paramin_l(1)>theta_l(t)||theta_l(t)>paramax_l(1)
        n = floor((theta_l(t)-paramin_l(1))/range(1));
        if mod(n,2)==1
            theta_l(t)=theta_l(t)-n*range(1);
        else
            theta_l(t)=paramax_l(1)+n*range(1)-(theta_l(t)-paramin_l(1));
        end
    end
    if paramin_l(2)>c_l(t)||c_l(t)>paramax_l(2)
        n = floor((c_l(t)-paramin_l(2))/range(2));
        if mod(n,2)==1
            c_l(t)=c_l(t)-n*range(2);
        else
            c_l(t)=paramax_l(2)+n*range(2)-(c_l(t)-paramin_l(2));
        end
    end
    if paramin_l(3)>g_l(t)||g_l(t)>paramax_l(3)
        n = floor((g_l(t)-paramin_l(3))/range(3));
        if mod(n,2)==1
            g_l(t)=g_l(t)-n*range(3);
        else
            g_l(t)=paramax_l(3)+n*range(3)-(g_l(t)-paramin_l(3));
        end
    end
    if paramin_l(4)>c_t(t)||c_t(t)>paramax_l(4)
        n = floor((c_t(t)-paramin_l(4))/range(4));
        if mod(n,2)==1
            c_t(t)=c_t(t)-n*range(4);
        else
            c_t(t)=paramax_l(4)+n*range(4)-(c_t(m)-paramin_l(4));
        end
    end
    if paramin_l(5)>g_t(t)||g_t(t)>paramax_l(5)
        n = floor((g_t(t)-paramin_l(5))/range(5));
        if mod(n,2)==1
            g_t(t)=g_t(t)-n*range(5);
        else
            g_t(m)=paramax_l(5)+n*range(5)-(g_t(m)-paramin_l(5));
        end
    end
    
    % Default Intensity Model Results Visualization
    options = optimset('Display','off');
    lambda0(t) = fsolve(@(x)CDXPIDE(0,0,r0(t),abs(x),[dT:dT:CDXMat],...
         CDX_0*1e-4,delta,U,theta_r(t),theta_l(t),rho(t),c_r(t),c_l(t),c_t(t),...
         g_r(t),g_l(t),g_t(t),0,0,0),0,options);
    fprintf('lambda0 = %d; theta_l = %d; rho = %d; c_l = %d; g_l = %d; c_t = %d; g_t = %d;\n',lambda0(t),theta_l(t),rho(t),c_l(t),g_l(t),c_t(t),g_t(t))
    fprintf('CDX Market calibration result\n')
    CDXOMCCalibRecPay(r0(t),delta,U,Gr,seed,...
        theta_r(t),theta_l(t),rho(t),c_r(t),c_l(t),c_t(t),g_r(t),g_l(t),g_t(t),Data_t,T0,ind,dT,CDXMat,1,CDX_0,paramin_l,paramax_l,1,...
        0); 
        
    %Market Implied Spread Moments
    payer_t_mkt = MidPay*1e+4;
    receiver_t_mkt = MidRec*1e+4;
    strike_t = kappa*1e+4;
    c_hat(t) = strike_t(end);
    c_f_mkt(t) = fsolve(@(x)interp1(strike_t,payer_t_mkt-receiver_t_mkt,x),CDX_0,options);
    A_0_mkt(t) = interp1(strike_t,payer_t_mkt-receiver_t_mkt,0,'linear','extrap')/c_f_mkt(t);
    cmin = min(strike_t);
    cmax = max(strike_t);
    
    [~,kappa_f_ind]=min((strike_t-c_f_mkt(t)).^2); %Assume chat = c_f
    if strike_t(kappa_f_ind)>c_f_mkt(t)
        strike_t_receiver = [strike_t(1:kappa_f_ind-1);c_f_mkt(t)];
        strike_t_payer = [c_f_mkt(t);strike_t(kappa_f_ind:end)];
    else
        strike_t_receiver = [strike_t(1:kappa_f_ind);c_f_mkt(t)];
        strike_t_payer = [c_f_mkt(t);strike_t(kappa_f_ind+1:end)];
    end
    
    receiver_t_mkt = interp1(strike_t,receiver_t_mkt,strike_t_receiver,'linear','extrap');
    payer_t_mkt = interp1(strike_t,payer_t_mkt,strike_t_payer,'linear','extrap');
    
    Vol(t) = trapz(strike_t_receiver,2*receiver_t_mkt)...
              +trapz(strike_t_payer,2*payer_t_mkt);
    Cub(t) = trapz(strike_t_receiver,6*(strike_t_receiver-c_f_mkt(t)).*receiver_t_mkt)...
              +trapz(strike_t_payer,6*(strike_t_payer-c_f_mkt(t)).*payer_t_mkt);
    Qua(t) = trapz(strike_t_receiver,(12*(strike_t_receiver-c_f_mkt(t)).^2).*receiver_t_mkt)...
              +trapz(strike_t_payer,(12*(strike_t_payer-c_f_mkt(t)).^2).*payer_t_mkt);
        
	Var(t) = Vol(t)/A_0_mkt(t);
    Skew(t) = ( Cub(t)/A_0_mkt(t) ) / ( Var(t)^(3/2) );
    Kurt(t) = ( Qua(t)/A_0_mkt(t) ) / ( Var(t)^(2) );
    
    %Model Implied Option Surface 
    rng('default'); rng(seed);
    Gt = gamrnd(g_t(t)*dt,1/c_t(t),Nsim,Nsim);
    Gl = gamrnd(g_l(t)*Gt,1/c_l(t),Nsim,Nsim);
    er = exp( -theta_r(t) * ( T0-cumsum(ones(Nsim,1))*dt ) );
    el = exp( -theta_l(t) * ( T0-cumsum(ones(Nsim,1))*dt ) );
    eY_r = ( 1 - exp( -theta_r(t) * ( T0-cumsum(ones(Nsim,1))*dt ) ) ) / theta_r(t);
    rT0 = r0(t) * exp ( -theta_r(t)*T0 ) + Gr*er; %simulated short rate @ T0(m)
    lT0 = lambda0(t) * exp( -theta_l(t)*T0 ) + (rho(t)*Gr+Gl)*el; %simulated default intensity @ T0(m)
    Y_rT0 = r0(t) * ( 1-exp(-theta_r(t)*T0) ) / theta_r(t) + Gr*eY_r; %simulated integrated short rate @ T0(m)
    
    strike_t_mdl = kappa;
    payer_t_mdl = zeros(length(strike_t_mdl),1);
    receiver_t_mdl = zeros(length(strike_t_mdl),1);
    for kk = 1:length(strike_t_mdl)
        for s=1:Nsim
            payer_t_mdl(kk) = payer_t_mdl(kk) +  exp(-Y_rT0(s))...
                *max(-CDXPIDE(T0,T0,rT0(s),lT0(s),[T0+dT:dT:CDXMat],strike_t_mdl(kk),delta,U,...
                theta_r(t),theta_l(t),rho(t),c_r(t),c_l(t),c_t(t),g_r(t),g_l(t),g_t(t),0,0,0),0)/Nsim;
            receiver_t_mdl(kk) = receiver_t_mdl(kk) +  exp(-Y_rT0(s))...
                *max(CDXPIDE(T0,T0,rT0(s),lT0(s),[T0+dT:dT:CDXMat],strike_t_mdl(kk),delta,U,...
                theta_r(t),theta_l(t),rho(t),c_r(t),c_l(t),c_t(t),g_r(t),g_l(t),g_t(t),0,0,0),0)/Nsim;
        end
    end
    
    % Model Implied Annuity Price
    payer_t_mdl = payer_t_mdl*1e+4;
    receiver_t_mdl = receiver_t_mdl*1e+4;
    strike_t_mdl = strike_t_mdl*1e+4;
    c_f_mdl(t) = fsolve(@(x)interp1(strike_t_mdl,payer_t_mdl-receiver_t_mdl,x,'linear','extrap'),CDX_0,options);
    A_0_mdl(t) = interp1(strike_t_mdl,payer_t_mdl-receiver_t_mdl,0,'linear','extrap')/c_f_mdl(t);
    
    % Model Implied Spread Moments
    
    receiver_t_mdl = interp1(strike_t_mdl,receiver_t_mdl,strike_t_receiver,'linear','extrap');
    payer_t_mdl = interp1(strike_t_mdl,payer_t_mdl,strike_t_payer,'linear','extrap');
    
    Vol_Model(t) = trapz(strike_t_receiver,2*receiver_t_mdl)...
              +trapz(strike_t_payer,2*payer_t_mdl);
    Cub_Model(t) = trapz(strike_t_receiver,6*(strike_t_receiver-c_f_mdl(t)).*receiver_t_mdl)...
              +trapz(strike_t_payer,6*(strike_t_payer-c_f_mdl(t)).*payer_t_mdl);
    Qua_Model(t) = trapz(strike_t_receiver,(12*(strike_t_receiver-c_f_mdl(t)).^2).*receiver_t_mdl)...
              +trapz(strike_t_payer,(12*(strike_t_payer-c_f_mdl(t)).^2).*payer_t_mdl);
    
    Var_Model(t) = Vol_Model(t) / A_0_mdl(t);
    Skew_Model(t) = ( Cub_Model(t)/A_0_mdl(t) ) / (Var_Model(t)^(3/2));
    Kurt_Model(t) = ( Qua_Model(t) / A_0_mdl(t) ) / (Var_Model(t)^(2));
    
    % Print results
      
    fprintf('\n Var, Skew, Kurt comparison\n')
    fprintf('Variance Mkt = %d, Variance Model = %d\n', Var(t), Var_Model(t));
    fprintf('Skewness Mkt = %d, Skewness Model = %d\n', Skew(t), Skew_Model(t));
    fprintf('kurtosis Mkt = %d, Kurtosis Model = %d\n', Kurt(t), Kurt_Model(t));
end

%% Variables Definition
% t = length(t0_l);
t0 = datetime(1900+year(t0_l(1:t)-1),month(t0_l(1:t)-1),day(t0_l(1:t)-1));

Vol_Model = Vol_Model(1:t)';
Vol = Vol(1:t)';
Cub_Model = Cub_Model(1:t)';
Cub = Cub(1:t)';
Qua_Model = Qua_Model(1:t)';
Qua = Qua(1:t)';

c_f_mdl = c_f_mdl(1:t)';
c_f_mkt = c_f_mkt(1:t)';
A_0_mdl = A_0_mdl(1:t)';
A_0_mkt = A_0_mkt(1:t)';

c_hat = c_hat(1:t)';

Var_Model = Var_Model(1:t)';
Var = Var(1:t)';
Skew_Model = Skew_Model(1:t)';
Skew = Skew(1:t)';
Kurt_Model = Kurt_Model(1:t)';
Kurt = Kurt(1:t)';

% Parameters
r0 = r0(1:t)';
theta_r = theta_r(1:t)';
c_r = c_r(1:t)';
g_r = g_r(1:t)';
%%
lambda0 = [0;lambda0(1:t-1)'];
theta_l = [3.353399e+00;theta_l(1:t-1)'];
rho = [1.548835e-01;rho(1:t-1)'];
c_l = [4.317856e+00;c_l(1:t-1)'];
g_l = [6.061738e+00;g_l(1:t-1)'];
c_t = [1.900001e+02;c_t(1:t-1)'];
g_t = [3.529832e+00;g_t(1:t-1)'];

%% Save variables
A_0 = [A_0_mkt,A_0_mdl];
c_f = [c_f_mkt,c_f_mdl];
param = [r0,theta_r,c_r,g_r,lambda0,theta_l,c_l,g_l,c_t,g_t,rho];
VolStat = [Vol,Vol_Model];
CubStat = [Cub,Cub_Model];
QuaStat = [Qua,Qua_Model];

t0 = t0(1:end);
c_hat = c_hat(1:end);

ssave = 1;

if ssave == 1
    save(strcat('Annuity',num2str(CalibrationMaturity),'.mat'),'A_0')
    save(strcat('forwardspread',num2str(CalibrationMaturity),'.mat'),'c_f')
    save(strcat('Param',num2str(CalibrationMaturity),'.mat'),'param')
    save(strcat('VolContract',num2str(CalibrationMaturity),'.mat'),'VolStat')
    save(strcat('CubContract',num2str(CalibrationMaturity),'.mat'),'CubStat')
    save(strcat('QuaContract',num2str(CalibrationMaturity),'.mat'),'QuaStat')
    save(strcat('today',num2str(CalibrationMaturity),'.mat'),'t0')
    save(strcat('maxstrike',num2str(CalibrationMaturity),'.mat'),'c_hat')
end

%% Load variables

lload = 1;

if lload == 1
    clear
    clc
    CalibrationMaturity = 1;
    load(strcat('Annuity',num2str(CalibrationMaturity),'.mat'),'A_0')
    load(strcat('forwardspread',num2str(CalibrationMaturity),'.mat'),'c_f')
    load(strcat('Param',num2str(CalibrationMaturity),'.mat'),'param')
    load(strcat('VolContract',num2str(CalibrationMaturity),'.mat'),'VolStat')
    load(strcat('CubContract',num2str(CalibrationMaturity),'.mat'),'CubStat')
    load(strcat('QuaContract',num2str(CalibrationMaturity),'.mat'),'QuaStat')
    load(strcat('today',num2str(CalibrationMaturity),'.mat'),'t0')
    load(strcat('maxstrike',num2str(CalibrationMaturity),'.mat'),'c_hat')
    
    Vol=VolStat(:,1);
    Vol_Model=VolStat(:,2);
    Cub=CubStat(:,1);
    Cub_Model=CubStat(:,2);
    Qua=QuaStat(:,1);
    Qua_Model=QuaStat(:,2);
    
    r0 = param(:,1);
    lambda0 = param(:,5);
    rho = param(:,11);
    
    c_f_mkt = c_f(:,1);
    c_f_mdl = c_f(:,2);
    
    A_0_mkt = A_0(:,1);
    A_0_mdl = A_0(:,2);
    AA_0 = A_0;
    
    VarStat=VolStat./AA_0;
    Var=VarStat(:,1);
    Var_Model=VarStat(:,2);
    SkewStat=(CubStat./AA_0)./(VarStat.^(3/2));
    Skew=SkewStat(:,1);
    Skew_Model=SkewStat(:,2);
    KurtStat=(QuaStat./AA_0)./(VarStat.^(2));
    Kurt=KurtStat(:,1);
    Kurt_Model=KurtStat(:,2);

    close all
end

%% Results Var, Skew, Kurt

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
box on
hold on
grid on
plot(t0(2:end),Var(2:end))
plot(t0(2:end),Var_Model(2:end))
set(gca,'YScale','log')
legend('Market','Model','Interpreter','latex','Location','Southeast')
title('Variance')
hold off
str=strcat('Var',num2str(CalibrationMaturity));
fname=str;
saveas(gcf, fullfile(fpath, fname), 'epsc');

figure
box on
hold on
grid on
plot(t0(2:end),Skew(2:end))
plot(t0(2:end),Skew_Model(2:end))
set(gca,'YScale','log')
legend('Market','Model','Interpreter','latex','Location','Southeast')
title('Skewness')
hold off
str=strcat('Skew',num2str(CalibrationMaturity));
fname=str;
saveas(gcf, fullfile(fpath, fname), 'epsc');

figure
box on
hold on
grid on
plot(t0(2:end),Kurt(2:end))
plot(t0(2:end),Kurt_Model(2:end))
set(gca,'YScale','log')
legend('Market','Model','Interpreter','latex','Location','Southeast')
title('Kurtosis')
hold off
str=strcat('Kur',num2str(CalibrationMaturity));
fname=str;
saveas(gcf, fullfile(fpath, fname), 'epsc');

rhovar = corrcoef(Var(2:end-1),Var_Model(2:end-1));
rhoskew = corrcoef(Skew(2:end-1),Skew_Model(2:end-1));
rhokurt = corrcoef(Kurt(2:end-1),Kurt_Model(2:end-1));

% Print
fprintf('\nCorrelation Results\n')
fprintf('Correlation Implied/Model Spread Variance: %d\n', rhovar(1,2))
fprintf('Correlation Implied/Model Spread Skew: %d\n', rhoskew(1,2))
fprintf('Correlation Implied/Model Spread kurtosis: %d\n', rhokurt(1,2))

%% Results: daily short rate, default intensity and correlation
figure
hold on
box on
grid on
plot(t0(2:end),r0(2:end))
ylabel('$r(t)$','interpreter','latex')
hold off
str=strcat('r',num2str(CalibrationMaturity));
fname=str;
saveas(gcf, fullfile(fpath, fname), 'epsc');

figure
hold on
box on
grid on
plot(t0(2:end),lambda0(2:end))
ylabel('$\lambda(t)$','interpreter','latex')
hold off
str=strcat('lambda',num2str(CalibrationMaturity));
fname=str;
saveas(gcf, fullfile(fpath, fname), 'epsc');

figure
hold on
box on
grid on
plot(t0(2:end),rho(2:end))
ylabel('$\rho(t)$','interpreter','latex')
hold off
str=strcat('rho',num2str(CalibrationMaturity));
fname=str;
saveas(gcf, fullfile(fpath, fname), 'epsc');