%% Fourier Inversion: Short Rate
clear
clc

% Model Paramaters
t = 1;

% Model Paramaters
r0 = 0.01468981; theta_r = [0.1,0.5500511,2.5]; c_r = 400.0005; g_r = 3.947560;

Nsim = 100;
rng('default'); %freeze the seed
rng(1);
U = rand(Nsim,1);

B = 1000000;
N=2^20;
eta=B/N;
lambda = 2*pi/B;
b = 0;%lambda*N/2;
alpha=[0:N-1]*eta;
w = ones(1,N); w(1)=1/2;
r = -b+lambda*[0:N-1];

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
for j=1:length(theta_r)
    phi = exp( -g_r*sum( log(1-(1i/c_r)*alpha.*exp(-theta_r(j)*(t-t*U))) )*t/Nsim );
    f = fft((1/pi)*exp(1i*alpha*b).*phi*eta);
    [~,Nmax] = min((r-2*r0).^2);
    %Nmax=N;
    plot(r(1:Nmax),f(1:Nmax));
    str{j} = strcat('$\theta_r = $',num2str(theta_r(j)));
end
legend (str,'interpreter','latex');
hold off
str=strcat('r_Density');
fname=str;
saveas(gcf, fullfile(fpath, fname), 'epsc');
close all