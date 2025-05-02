%% FFT Gamma density

clear
clc

alpha = 2;
beta = 0.5;
mu = 1;
sigma = 1;

B = 1000;
N=2^15;
eta=B/N;
lambda = 2*pi/B;
b = 0;%lambda*N/2;
u=[0:N-1]*eta;
w = ones(1,N); w(1)=1/2;
x = -b+lambda*[0:N-1];
phi = ((1-1i*u/beta).^(-alpha)).*w;%exp(1i*mu*u-(sigma*u).^2/2).*w;%
f = fft((1/pi)*exp(1i*u*b).*phi*eta);
f_ex=((beta^alpha)/gamma(alpha))*x.^(alpha-1).*exp(-beta*x);%(1/(sigma*sqrt(2*pi)))*exp(-0.5*((x-mu)/sigma).^2);%

figure
hold on
grid on
box on
[~,Nmax] = min((x-50).^2);
%Nmax=N;
plot(x(1:Nmax),f(1:Nmax));
plot(x(1:Nmax),f_ex(1:Nmax))
legend('Numerical','Analytical')
hold off