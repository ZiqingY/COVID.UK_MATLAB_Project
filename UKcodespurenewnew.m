set(groot,'defaultLineLineWidth',2.0);
zeta = 1/30;
%zeta = 10;
gamma = 0.4;
sigma = 0.425;
eta = 0.025;
nu = 0.2;
betabar = 3.6*gamma;
betabarv = 5*gamma;
kappa = 1.4000e+05;

population = 67900000;

seasonalsize = 1;
seasonalposition = 15;
fatiguesize = 0.575;
fatiguemean = 305;
fatiguesig = 8;

tv = 191;
Evbar = 1/population;

E0 = 10/population;
Ev0 = 0;
S0 = 1-E0;
I0 = 0;
Iv0 = 0;
R0 = 0;
H0 = 0;
D0 = 0;

cases = 1;

t0 = 0;  
tfinal = 2*356;

y0 = zeros(23,1);
y0(1) = S0;
y0(2) = E0;
y0(3) = Ev0;
y0(4) = I0;
y0(5) = Iv0;
y0(6) = R0;
y0(7) = H0;
y0(8) = D0;
y0(9) = betabar;
y0(10) = kappa;
y0(11) = zeta;
y0(12) = cases;
y0(13) = gamma;
y0(14) = sigma;
y0(15) = eta;
y0(16) = nu;
y0(17) = seasonalsize;
y0(18) = fatiguesize;
y0(19) = fatiguemean;
y0(20) = fatiguesig;
y0(21) = betabarv;
y0(22) = tv;
y0(23) = Evbar;
y0(24) = seasonalposition;

opts = odeset('Reltol',1e-13,'AbsTol',1e-14,'Stats','on');
[t1,y1] = ode113(@bsirSEIRHD,[t0 tfinal],y0,opts);


S1 = y1(:,1);
E1 = y1(:,2);
Ev1 = y1(:,3);
I1 = y1(:,4);
Iv1 = y1(:,5);
R1 = y1(:,6);
H1 = y1(:,7);
D1 = y1(:,8);
Ddot1 = nu*zeta*(eta*gamma*(I1+Iv1) - zeta*H1);
gD1 = Ddot1./(nu.*zeta.*H1);

load ukdatanewnew.mat
tdata = ukdatanewnew(:,1)-ukdatanewnew(1,1);
t1date = datetime('15-Feb-2020')+t1;
tdatadate = datetime('15-Feb-2020')+tdata;

figure(1)
% plot death
subplot(2,1,1)
plot(t1date,nu*zeta*H1*population,tdatadate,ukdatanewnew(:,4))
line([37 37],[0 1500],'linestyle', '--' , 'Color','r', 'LineWidth', 1)
line([129 129],[0 1500],'linestyle', '--' , 'Color','g', 'LineWidth', 1)
line([259 259],[0 1500],'linestyle', '--' , 'Color','r', 'LineWidth', 1)
line([291 291],[0 1500],'linestyle', '--' , 'Color','g', 'LineWidth', 1)
line([324 324],[0 1500],'linestyle', '--' , 'Color','r', 'LineWidth', 1)
line([422 422],[0 1500],'linestyle', '--' , 'Color','g', 'LineWidth', 1)

% plot fatigue
subplot(2,1,2)
tnew=transpose((1:570))
knew = kappa*(1 -normcdf(tnew,fatiguemean,fatiguesig))+fatiguesize*kappa*normcdf(tnew,fatiguemean,fatiguesig);
t1datenew = datetime('15-Feb-2020')+tnew
plot(t1datenew, knew)
line([37 37],[70000 150000],'linestyle', '--' , 'Color','r', 'LineWidth', 1)
line([129 129],[70000 150000],'linestyle', '--' , 'Color','g', 'LineWidth', 1)
line([259 259],[70000 150000],'linestyle', '--' , 'Color','r', 'LineWidth', 1)
line([291 291],[70000 150000],'linestyle', '--' , 'Color','g', 'LineWidth', 1)
line([324 324],[70000 150000],'linestyle', '--' , 'Color','r', 'LineWidth', 1)
line([422 422],[70000 150000],'linestyle', '--' , 'Color','g', 'LineWidth', 1)

figure(2)
% plot hospitalization
plot(t1date,nu*zeta*H1*population,tdatadate,ukdatanewnew(:,4))
