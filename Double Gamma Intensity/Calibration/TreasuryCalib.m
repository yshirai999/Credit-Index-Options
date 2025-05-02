function P = TreasuryCalib(r,theta_r,c_r,g_r,T,Y,U,w,min,max,ellmax,Display)
range = max-min;
if min(1)>theta_r||theta_r>max(1)
    n = floor((theta_r-min(1))/range(1));
    if mod(n,2)==1
        theta_r=theta_r-n*range(1);
    else
        theta_r=max(1)+n*range(1)-(theta_r-min(1));
    end
end

if min(2)>c_r||c_r>max(2)
    n = floor((c_r-min(2))/range(2));
    if mod(n,2)==1
        c_r=c_r-n*range(2);
    else
        c_r=max(2)+n*range(2)-(c_r-min(2));
    end
end

if min(3)>g_r||g_r>max(3)
    n = floor((g_r-min(3))/range(3));
    if mod(n,2)==1
        g_r=g_r-n*range(3);
    else
        g_r=max(3)+n*range(3)-(g_r-min(3));
    end
end

L = length(T);
DCF = zeros(1,L);
P = 0;
for ell = 1:L
    if T(ell)<=ellmax
        N = floor(T(ell)/0.5);
        for n=1:N
            DCF(ell) = DCF(ell) + ZCB(r,theta_r,g_r,c_r,0,0.5*n,U)*((1+Y(ell)/200)^(2*0.5)-1);
        end
        DCF(ell) = DCF(ell) + ZCB(r,theta_r,g_r,c_r,0,T(ell),U)*(1+Y(ell)/200)^(2*(T(ell)-0.5*N));
        if Display==1
            fprintf('Maturity = %4.2f, Model Price = %4.2f\n',T(ell),DCF(ell));
        end
        P = P + w(ell)*(DCF(ell)-1)^2;
    end
end

end

