function w = GTCInt(y,c,g,Y)
% Gamma Levy density weight for gamma time changed gamma process
Nsim = length(Y);
w =  sum( ( ((c*y').^(g*Y))./gamma(g*Y) ) ./ Y, 2 ) / Nsim;
end