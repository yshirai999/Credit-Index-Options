%% FFT Bivariate Stationary Density

clear
clc

%% Parameters

% Model Paramaters
theta_r = 0.5500511; c_r = 400.0005; g_r = 3.947560;
lambda0 = 0; theta_l = 7.869028e-01; rho = 1.562220e-01;
    c_l = 2.032924e+01; g_l = 4.122366e+00; c_t = 6.040000e+02; g_t = 3.319264e+00;

% Montecarlo integration parameters
Nsim = 100;
rng('default'); %freeze the seed
rng(1);
U = rand(Nsim,1);

%% Fourier Inversion

% Parameters
B = 750000;
N=2^13;
eta=2*B/N;
lambda = pi/B;
b = pi/eta;

% Frequency Grid
u1 = repmat(-B+[0:N-1]*eta,1,N);
u2 = kron(-B+[0:N-1]*eta,ones(1,N));
u=[u1;u2];

% Space Grid 
x1 = repmat(-b+[0:N-1]*lambda,1,N);
x2 = kron(-b+[0:N-1]*lambda,ones(1,N));
x=[x1;x2];

% Fourier Transform
int_r = zeros(1,N*N);
int_l = zeros(1,N*N);
t = 100;

for n=1:Nsim
    int_r = int_r + log ( 1 + (1i/c_r)...
        * ( exp(-theta_r*(t-t*U(n)))*u1 + rho*exp(-theta_l*(t-t*U(n)))*u2 )...
        ) * t/Nsim;
    int_l = int_l + log ( 1 + (g_l/c_t) * log( 1 + (1i/c_l) * exp(-theta_l*(t-t*U(n)) * u2 ) ) ) * t/Nsim;
    n
end
%     int_r = sum(log ( 1 + (1i/c_r)...
%         * ( exp(-theta_r*(t-t*U))*u1 + rho*exp(-theta_l*(t-t*U))*u2 )...
%         ) )* t/Nsim;
%     int_l = sum( log ( 1 + (1i/c_l) * exp(-theta_l*(t-t*U) * u2 ) ) )* t/Nsim;
phi = reshape( exp ( - g_r*int_r - g_t*int_l ), N, N );
H = (-1).^([0:N-1]+[0:N-1]').*phi;
f = ( (eta*N/(2*pi))^2) * (-1).^([0:N-1]+[0:N-1]') .* ifft2(H)';

%% Visualization

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

[~,Nmin_l] = min((x(1,1:N)-0.00025).^2);
[~,Nmax_l] = min((x(1,1:N)-0.0018).^2);
[~,Nmin_r] = min((x(1,1:N)-0.005).^2);
[~,Nmax_r] = min((x(1,1:N)-0.018).^2);
[xx,yy]=meshgrid(x(1,Nmin_r:Nmax_r),x(1,Nmin_l:Nmax_l));
zz=max(real(f(Nmin_l:Nmax_l,Nmin_r:Nmax_r)),0);

% [xx,yy]=meshgrid(x1(1:N),x1(1:N));
% zz=max(real(f),0);

figure(1)
hold on
grid on
box on
waterfall(xx,yy,zz)
xlabel('r')
ylabel('$\lambda$','interpreter','latex')
view(100,30)
hold off
str=strcat('BivDensityStationary');
fname=str;
saveas(gcf, fullfile(fpath, fname), 'epsc');
%close all

figure(2)
hold on
contour(xx,yy,zz)
xlabel('r')
ylabel('$\lambda$','interpreter','latex')
colorbar
hold off
str=strcat('BivDensityContourStationary');
fname=str;
saveas(gcf, fullfile(fpath, fname), 'epsc');
%close all