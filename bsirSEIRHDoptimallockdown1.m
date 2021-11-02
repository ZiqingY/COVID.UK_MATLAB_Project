function yp = bsirSEIRHDassumptionset2(t,y,params);

S = y(1);
E = y(2);
Ev = y(3);
I = y(4);
Iv = y(5);
R = y(6);
H = y(7);
D = y(8);
betabar= params(1);
kappa = params(2)
kappa1 = params(3);
kappa2 = params(4);
kappa3 = params(5);
zeta = params(6);
cases = params(7);
gamma = params(8);
sigma = params(9);
eta = params(10);
nu = params(11);
seasonalsize = params(12);
fatiguesize1 = params(13);
fatiguemean1 = params(14);
fatiguesig1 = params(15);
fatiguesize2 = params(16);
fatiguemean2 = params(17);
fatiguesig2 = params(18);
fatiguesize3 = params(19);
fatiguemean3 = params(20);
fatiguesig3 = params(21);
betabarv = params(22);
tv = params(23);
Evbar = params(24);
seasonalposition = params(25);


if cases == 1;
    psi = seasonalsize*(cos((t+seasonalposition)*2*pi/365)-1)/2;
    kt = (1 + (fatiguesize1-1)*normcdf(t,fatiguemean1,fatiguesig1))*kappa1.*(t<=102) + (1 + (fatiguesize2-1)*normcdf(t,fatiguemean2,fatiguesig2))*kappa2.*(t>102 & t<=182) + (1 + (fatiguesize3-1)*normcdf(t,fatiguemean3,fatiguesig3))*kappa3.*(t>182)
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


beta = betabar*exp(-kt*nu*zeta*H+psi);
betav = betabarv*exp(-kt*nu*zeta*H+psi);


x(1) = -beta*S*I - betav*S*Iv;
x(2) = beta*S*I - sigma*E;
x(3) = betav*S*Iv - sigma*Ev +(t<(tv+2))*(t>tv)*Evbar;
x(4) = sigma*E - gamma*I;
x(5) = sigma*Ev - gamma*Iv;
x(6) = (1-nu)*zeta*H + (1-eta)*gamma*(I+Iv) -(t<(tv+2))*(t>tv)*Evbar;
x(7) = eta*gamma*(I+Iv) - zeta*H;
x(8) = nu*zeta*H;

yp = transpose(x);