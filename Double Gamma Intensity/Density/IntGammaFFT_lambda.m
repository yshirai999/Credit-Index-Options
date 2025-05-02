%% Fourier Inversion: Default Intensity
clear
clc

%% Model Paramaters
t = 1;

% Model Paramaters
r0 = 0.01468981; theta_r = 0.5500511; c_r = 400.0005; g_r = 3.947560;
lambda0 = 0; theta_l = [0.1,7.869028e-01,4]; rho = 1.562220e-01;
    c_l = 2.032924e+01; g_l = 4.122366e+00; c_t = 6.040000e+02; g_t = 3.319264e+00;
    
Nsim = 100;
rng('default'); %freeze the seed
rng(1);
U = rand(Nsim,1);

%% Fourier Inversion
B = 1000000;
N=2^20;
eta=B/N;
lambda = 2*pi/B;
b = lambda*N/2;
beta=[0:N-1]*eta;
w = ones(1,N); w(1)=1/2;
l = -b+lambda*[0:N-1];

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
hold on
box on
grid on
phi = zeros(1,2*(N-1)+1);
[~,Nmax] = min((l-0.006).^2);
[~,Nmin] = min((l-0).^2);
for j=1:length(theta_l)
    phi = exp( 1i*beta*lambda0*exp(-theta_l(j)*t)...
        -g_r*sum( log(1-(1i/c_r)*rho*beta.*exp(-theta_l(j)*(t-t*U))) )*t/Nsim...
        -g_t*sum( log(1+(g_l/c_t)*log(1-(1i/c_l)*beta.*exp(-theta_l(j)*(t-t*U))) ))*t/Nsim);
    f = fft((1/pi)*exp(1i*beta*b).*phi*eta);
    plot(l(Nmin:Nmax),max(real(f(Nmin:Nmax)),0));
    str{j} = strcat('$\theta_{\lambda} = $',num2str(theta_l(j)));
end
legend (str,'interpreter','latex');
ax = gca;
ax.XAxis.Exponent = 0;
hold off
str=strcat('l_Density');
fname=str;
saveas(gcf, fullfile(fpath, fname), 'epsc');
close all