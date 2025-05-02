%% CDXO Price
clear
%clc
close all
%% Parameters
% CDX Parameters
T0 = [0.036111111111111, 0.133333333333333, 0.211111111111111, 0.288888888888889, 0.386111111111111, 0.463888888888889];
T0 = T0(2);
dT = 0.5;
T = T0+dT:dT:5;
kappa = 0.006;
delta = 0.6;

% Model Paramaters
r0 = 1.468981e-02; theta_r = 5.500511e-01; c_r = 4.000005e+02; g_r = 3.947560e+00;
lambda0 = 0; theta_l = 3.353399e+00; rho = 1.548835e-01;
    c_l = 4.317856e+00; g_l = 6.061738e+00; c_t = 1.900001e+02; g_t = 3.529832e+00;

% Grid Paramaters
rmin = 0;
rmax = 0.35;
lmin = 0;
lmax = rho*rmax;
M = 100;
dt = T0/M;
t = 0:dt:T0;

N = [50,100,150,200,250];
errMC = zeros(1,length(N));

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

for nn = 1:length(N)
    tic
    dr = rmax/N(nn);
    r = rmin:dr:rmax;
    
    L = N(nn);
    dl = lmax/L;
    lambda = lmin:dl:lmax;
    
    u = zeros((N(nn)-1)*(L-1),M+1);
    
    tic 
    % Initial condition (Payoff)
    Nsim = 100; rng('default'); rng(1); U = rand(1,Nsim);
    Y_r = sort(-log(U)/c_r); %Exponential variables to calculate integrals in d\phi_r
    Y_l = sort(-log(U)/c_l); %Exponential variables to calculate integrals in d\phi_l
    Y_t = sort(-log(U)/c_t); %Exponential variables to calculate integrals in d\phi_t
    w_tau = GTCInt(Y_t,c_l,g_l,Y_l); %weights for Gamma time changed gamma Levy measure
    u(:,1) = max(CDXPIDE(T0,T0,r(2:N(nn)),lambda(2:L),T,kappa,delta,U,...
          theta_r,theta_l,rho,c_r,c_l,c_t,g_r,g_l,g_t,N(nn)-2,L-2),0);
    
    %Evaluation of omega
    omega1 = 0;%(g_r/c_r);
    omega2 = 0;%rho*omega1+(g_l/c_l);
    
    % Diffusion Matrix
    % Block matrices
    c = (1+dt*(r(2:N(nn))))';
    w = (-dt*(-theta_r*r(3:N(nn)+1)-omega1)/(2*dr))';
    e = (dt*(-theta_r*r(1:N(nn)-1)-omega1)/(2*dr))';
    Ce = spdiags( [-w c -e], [-1, 0, 1], N(nn)-1, N(nn)-1 );
    
    % West Boundary condition
    w1 = (-dt*(-theta_r*r(2)-omega1)/(2*dr))';
    Ce(1,1)=Ce(1,1)-2*w1;
    Ce(1,2)=Ce(1,2)+w1;
    
    % East Boundary condition
    eN = (dt*(-theta_r*r(N(nn))-omega1)/(2*dr))';
    Ce(N(nn)-1,N(nn)-1)=Ce(N(nn)-1,N(nn)-1)-2*eN;
    Ce(N(nn)-1,N(nn)-2)=Ce(N(nn)-1,N(nn)-2)+eN;
    
    % Discrete Differential Operator
    ee = ones(N(nn)-1,1);
    I_L = speye(L-1);
    Lambda = spdiags([-lambda(3:L+1)'*theta_l-omega2,...
             lambda(1:L-1)'*theta_l+omega2],[-1 1],L-1,L-1)*dt/(2*dl);
         %+spdiags(lambda(2:L)',0,L-1,L-1)*dt;
    I_N = speye(N(nn)-1);
    A = kron(I_L,Ce)+kron(Lambda,I_N);
    
    % South Boundary condition
    for i=1:N(nn)-1
        A(i,i) = A(i,i) + 2 * dt*(-theta_l*lambda(2)-omega2)/(2*dl);
        A(i,N(nn)-1+i)=A(i,N(nn)-1+i) - dt*(-theta_l*lambda(2)-omega2)/(2*dl);
    end
    
    % North Boundary condition
    for i=1:N(nn)-1
         A((L-2)*(N(nn)-1)+i,(L-2)*(N(nn)-1)+i)...
            =A((L-2)*(N(nn)-1)+i,(L-2)*(N(nn)-1)+i)-2*dt*(-theta_l*lambda(L)-omega2)/(2*dl);
        A((L-2)*(N(nn)-1)+i,(L-3)*(N(nn)-1)+i)...
            =A((L-2)*(N(nn)-1)+i,(L-3)*(N(nn)-1)+i)+dt*(-theta_l*lambda(L)-omega2)/(2*dl);
    end
    
    % Finite Difference Scheme
    uu=zeros(L+1,N(nn)+1);
    uuinner=(reshape(u(:,1),L-1,N(nn)-1))';
    uu(2:L,2:N(nn)) = uuinner;
    uu(2:L,1) = 2*uu(2:L,2)-uu(2:L,3); %west boundary
    uu(2:L,N(nn)+1) = 2*uu(2:L,N(nn))-uu(2:L,N(nn)-1); %east boundary
    uu(L+1,1:N(nn)+1) = 2*uu(L,1:N(nn)+1)-uu(L-1,1:N(nn)+1); % north boundary
    uu(1,1:N(nn)+1) = 2*uu(2,1:N(nn)+1)-uu(3,1:N(nn)+1); %south boundary
    
    for j=2:M+1
        R=Rfun(r,lambda,uu,Y_r,Y_t,w_tau,Nsim,rho,g_r,g_t,c_r,c_t,N(nn),L);
        u(:,j)=A\(u(:,j-1)+dt*R);
        uuinner=(reshape(u(:,j),N(nn)-1,L-1))';
        uu(2:L,2:N(nn)) = uuinner;
        uu(2:L,1) = 2*uu(2:L,2)-uu(2:L,3); %west boundary
        uu(2:L,N(nn)+1) = 2*uu(2:L,N(nn))-uu(2:L,N(nn)-1); %east boundary
        uu(L+1,1:N(nn)+1) = 2*uu(L,1:N(nn)+1)-uu(L-1,1:N(nn)+1); %north boundary
        uu(1,1:N(nn)+1) = 2*uu(2,1:N(nn)+1)-uu(3,1:N(nn)+1); %south boundary
    end
    
%     CDXOprice1 = max(interp2(r,lambda,uu,r0,lambda1),0);
%     cputime = toc;
%     CDXOprice2 = max(interp2(r,lambda,uu,r0,lambda2),0);
%     CDXOprice3 = max(interp2(r,lambda,uu,r0,lambda3),0);
    
    CDXOprice0 = max(interp2(r,lambda,uu,r0,lambda0),0);
    
    cputime = toc;
    
    % Visualization
%     [~,nt] = min((lambda-lambda0).^2); %For higher values of lambda, one would need a bigger initial rmax, lambdamax.
%     range = 2^5;
%     nnt = [max(floor(nt/range),1):floor(nt*(range+1)/range)];
%     uu2 = uu(nnt,nnt);
%     rtest = r(nnt);
%     lambdatest = lambda(nnt);
%     [rr,ll]=meshgrid(rtest,lambdatest);
    
%     figure
%     box on
%     grid on
%     surf(rr,ll,uu2)
%     xlabel('$r$','Interpreter','Latex')
%     ylabel('$\lambda$','Interpreter','Latex')
%     view(160,45)
    
    Nsim = 100;
    dtsim = T0/Nsim;
    rng('default'); rng(1);
    Gr = gamrnd(g_r*dtsim,1/c_r,Nsim,Nsim);
    Gt = gamrnd(g_t*dtsim,1/c_t,Nsim,Nsim); %Gamma random time increments
    Gl = gamrnd(g_l*Gt,1/c_l,Nsim,Nsim);
    pMC=CDXOMC(0,T0,r,lambda,T,kappa,delta,U,Gr,Gl,...
        theta_r,theta_l,rho,c_r,c_l,c_t,g_r,g_l,g_t,N(nn),L);
    pMC=reshape(pMC,N(nn)+1,L+1)';
    
%     pMC2 = pMC(nnt);
%     err = pMC-uu;
    errMC(nn) = max(max(abs(uu-pMC)))*1e4;
    
    pMC0 = CDXOMC(0,T0,r0,lambda0,T,kappa,delta,U,Gr,Gl,...
        theta_r,theta_l,rho,c_r,c_l,c_t,g_r,g_l,g_t,0,0);
    
%     figure
%     box on
%     grid on
%     hold on
%     [~,nt_r] = min((r-r0).^2);
%     plot(lambda,pMC(:,nt_r))
%     plot(lambda,uu(:,nt_r))
%     xlabel('$\lambda$','Interpreter','Latex')
%     legend('Montecarlo','PIDE')
%     hold off

%      figure
%      box on
%      grid on
%      hold on
%      [rr,ll]=meshgrid(r,lambda);
%      surf(rr,ll,pMC)
%      xlabel('$r$','Interpreter','Latex')
%      ylabel('$\lambda$','Interpreter','Latex')
%      title('MC')
%      view(160,45)
%      hold off
     
    figure
    box on
    grid on
    hold on
    [rr,ll]=meshgrid(r,lambda);
    surf(rr,ll,uu*1e4)
    xlabel('$r$','Interpreter','Latex')
    ylabel('$\lambda$','Interpreter','Latex')
    view(160,45)
    hold off
    str=strcat('CDXO',num2str(N(nn)),num2str(M),num2str(Nsim));
    fname=str;
    saveas(gcf, fullfile(fpath, fname), 'epsc');
    %close all

    fprintf('N = %d\n',N(nn))
    fprintf('lambda0 = %d, CDXO price1 = %d, cputime = %d, pMC = %d\n',lambda0,CDXOprice0,cputime,pMC0)
%     fprintf('lambda1 = %d, CDXO price1 = %d, cputime = %d\n',lambda1,CDXOprice1,cputime)
%     fprintf('lambda2 = %d, CDXO price2 = %d, cputime = %d\n',lambda2,CDXOprice2,cputime)
%     fprintf('lambda3 = %d, CDXO price3 = %d, cputime = %d\n',lambda3,CDXOprice3,cputime)
end

% Visualization: MonteCarlo Error
figure
    
box on
grid on
plot(N,errMC,'linewidth',1.2,'marker','o')
xlabel('$N$','Interpreter','Latex')
set(gca,'YScale','log')
str=strcat('CDXOerrMC');
fname=str;
saveas(gcf, fullfile(fpath, fname), 'epsc');
close all

% Visualization: Spread
% figure(3)
% kappa = CDXSpread(0,T0,r,lambda,T,delta,U,...
%     theta_r,theta_l,rho,c_r,c_l,g_r,g_l,N(1),L);
% kappa=(reshape(kappa,N(1)+1,L+1))';
%     
% [rr,ll]=meshgrid(r,lambda);
% box on
% grid on
% surf(rr,ll,kappa)
% xlabel('$r$','Interpreter','Latex')
% ylabel('$\lambda$','Interpreter','Latex')
% view(45,45)
%     
% str=strcat('CDXSpread',num2str(N(1)),num2str(M),num2str(Nsim));
% fname=str;
% saveas(gcf, fullfile(fpath, fname), 'epsc');