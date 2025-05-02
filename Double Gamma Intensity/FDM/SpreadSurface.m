%% Spread Surface
clear
clc

rng('default'); %freeze the seed
rng(1);

%% Parameters
% CDX Parameters
T0 = 1;
dT = 0.25;
T = T0+dT:dT:5;
r0 = 0.0036;
delta = 0.6;

% Model Paramaters
theta_r = 0.18;
theta_l = 0.1;
rho = 0.25;
c_r = 571;
g_r = 4.98;
c_l = 653;
g_l = 7.04;
c_t = 4;
g_t = 400;

% Grid Paramaters
N = 100;
rmin = 0;
rmax = 0.15;
dr = rmax/N;
r = rmin:dr:rmax;

lmin = 0;
lmax = rho*rmax;
L = N;
dl = lmax/L;
lambda = lmin:dl:lmax;

M = 100;
dt = T0/M;
t = 0:dt:T0;

u = zeros((N-1)*(L-1),M+1);

tic 
% Spread surface (Payoff)
Nsim = 100;
U = rand(1,Nsim);
Y_r = -log(U)/c_r; %Exponential variables to calculate integrals in d\phi_r
Y_l = -log(U)/c_l; %Exponential variables to calculate integrals in d\phi_l

kappa = CDXSpread(0,T0,r,lambda,T,delta,U,...
    theta_r,theta_l,rho,c_r,c_l,c_t,g_r,g_l,g_t,N,L);
kappa=(reshape(kappa,N+1,L+1))';

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

figure(1)
[rr,ll]=meshgrid(r,lambda);
box on
grid on
surf(rr,ll,kappa)
xlabel('$r$','Interpreter','Latex')
ylabel('$\lambda$','Interpreter','Latex')
view(45,45)

str=strcat('CDXSpread',num2str(N),num2str(M),num2str(Nsim));
fname=str;
saveas(gcf, fullfile(fpath, fname), 'epsc');
close all


