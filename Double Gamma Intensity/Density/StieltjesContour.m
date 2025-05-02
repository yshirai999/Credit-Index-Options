%% Density
clear
clc

% Model Paramaters
t = 0.75;
theta_r = 0.18;
%c_r = 571;
g_r = 4.98;

% zeta function
zeta = @(s) exp( (1/theta_r)...
        * ( dilog( 1-(-1./s) ) - dilog( 1-(-exp(-theta_r*g_r*t)./s) )...
        - theta_r*g_r*t*log(s) ) );
    
zetaprime = @(s) zeta(s) .* exp( ( 1./(s*theta_r) )...
        .* ( log(1+1./s) - log(1+exp(-theta_r*g_r*t)./s)...
        - theta_r*g_r*t ) );

% integration
dy = 0.05;
y = dy:dy:1-dy;
q=zeros(length(y),1);
parfor k=1:length(y)
    fun = @(w) -(y(k)^t/(2*pi*1i))* ( (1+w).^(t-1)) .* zetaprime(y(k)*w );
    g = @(theta) cos(theta) + 1i*sin(theta);
    gprime = @(theta) -sin(theta) + 1i*cos(theta);
    q(k) = integral(@(theta) fun(g(theta)).*gprime(theta),-pi,pi);
end
%%
figure(1)
hold on
box on
grid on
plot(y,real(q)/(sum(real(q))*dy))
%legend ( 'density');
hold off
