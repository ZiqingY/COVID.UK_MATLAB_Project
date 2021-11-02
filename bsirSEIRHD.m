function yp = bsirSEIRHD(t,y);

S = y(1);
E = y(2);
Ev = y(3);
I = y(4);
Iv = y(5);
R = y(6);
H = y(7);
D = y(8);
betabar= y(9);
kappa = y(10);
zeta = y(11);
cases = y(12);
gamma = y(13);
sigma = y(14);
eta = y(15);
nu = y(16);
seasonalsize = y(17);
fatiguesize = y(18);
fatiguemean = y(19);
fatiguesig = y(20);
betabarv = y(21);
tv = y(22);
Evbar = y(23);
seasonalposition = y(24);


if cases == 1;
    psi = seasonalsize*(cos((t+seasonalposition)*2*pi/365)-1)/2;
    kappa = kappa*(1 -normcdf(t,fatiguemean,fatiguesig))+fatiguesize*kappa*normcdf(t,fatiguemean,fatiguesig);
end;

if cases == 2;
psi = 0;
end;

if cases == 3;
    psi = seasonalsize*(cos((t+seasonalposition)*2*pi/365)-1)/2;
end;

if cases == 4;
    psi = 0;
    kappa = kappa*(1 -normcdf(t,fatiguemean,fatiguesig))+fatiguesize*kappa*normcdf(t,fatiguemean,fatiguesig);
end;






beta = betabar*exp(-kappa*nu*zeta*H+psi);
betav = betabarv*exp(-kappa*nu*zeta*H+psi);



x(1) = -beta*S*I - betav*S*Iv;
x(2) = beta*S*I - sigma*E;
x(3) = betav*S*Iv - sigma*Ev +(t<(tv+2))*(t>tv)*Evbar;
x(4) = sigma*E - gamma*I;
x(5) = sigma*Ev - gamma*Iv;
x(6) = (1-nu)*zeta*H + (1-eta)*gamma*(I+Iv) -(t<(tv+2))*(t>tv)*Evbar;
x(7) = eta*gamma*(I+Iv) - zeta*H;
x(8) = nu*zeta*H;
x(9) = 0;
x(10) = 0;
x(11) = 0;
x(12) = 0;
x(13) = 0;
x(14) = 0;
x(15) = 0;
x(16) = 0;
x(17) = 0;
x(18) = 0;
x(19) = 0;
x(20) = 0;
x(21) = 0;
x(22) = 0;
x(23) = 0;
x(24) = 0;

yp = transpose(x);