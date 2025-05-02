%% Density
clear
clc

% Model Paramaters
t = 0.75;
theta_r = 0.18;
%c_r = 571;
g_r = 4.98;

dy = 0.05;
y = dy:dy:1-dy;
q=zeros(length(y),1);

Nsim = 500;
rng('default'); %freeze the seed
rng(1);
theta = -pi+2*pi*rand(1,Nsim);

parfor k=1:length(y)
    w = cos(theta) + 1i*sin(theta);
    zeta = exp( (1/theta_r)...
        * ( dilog( 1-(-1./(y(k)*w)) ) - dilog( 1-(-exp(-theta_r*g_r*t)./(y(k)*w))) )...
        - theta_r*g_r*t*log(y(k)*w) );
    
    zetaprime = zeta .* exp( ( 1./(y(k)*w*theta_r) )...
        .* ( log(1+1./(y(k)*w)) - log(1+exp(-theta_r*g_r*t)./(y(k)*w))...
        - theta_r*g_r*t ) );

    fun = -(y(k)^t/(2*pi*1i))* ( (1+w).^(t-1)) .* zetaprime;
    gprime = -sin(theta) + 1i*cos(theta);
    q(k) = sum(fun.*gprime)*(2*pi)/Nsim;
end
%%
figure(1)
hold on
box on
grid on
plot(y,real(q))
%legend ( 'density');
hold off