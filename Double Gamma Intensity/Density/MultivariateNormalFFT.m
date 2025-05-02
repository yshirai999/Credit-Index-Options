%% FFT Bivariate Normal density
clear
clc

%% Parameters
mu = [0.5,0];
sigma = [0.9,0.5;0.5,1.2];

%% Fourier Inversion

% Parameters
B = 500;
N=2^10;
eta=2*B/N;
lambda = pi/B;
b = lambda*N/2;%0;

% Frequency Grid
u1 = kron(-B+[0:N-1]*eta,ones(1,N));
u2 = repmat(-B+[0:N-1]*eta,1,N);
u=[u1;u2];

% Space Grid 
x1 = kron(-b+[0:N-1]*lambda,ones(1,N));
x2 = repmat(-b+[0:N-1]*lambda,1,N);
x=[x1;x2];

% Fourier Transform
phi = exp(1i*(mu*u)'-sum((u'*sigma).*u',2)/2);
H = eta*reshape(exp(1i*pi*sum(u))'.*phi,N,N)';
f = ( (eta*N/(2*pi))^2) * exp(1i*pi*([0:N-1]'+[0:N-1])) .* ifft2(H);

%% Exact Solution
f_ex = (1/(det(sigma)*2*pi))*exp( -0.5*( sum( ((x-mu')'/sigma) .* (x-mu')',2) ) );
f_ex = reshape(f_ex,N,N)';

%% Visualization
[xx,yy]=meshgrid(x(2,1:N),x(2,1:N));

figure
hold on
grid on
box on
waterfall(xx,yy,real(f(1:N,1:N)))
title('Numerical')
ylabel('r')
xlabel('$\lambda$','interpreter','latex')
view(45,45)
hold off

figure(2)
hold on
grid on
box on
waterfall(xx,yy,f_ex(1:N,1:N))
ylabel('r')
xlabel('$\lambda$','interpreter','latex')
title('Analytical')
view(45,45)
hold off