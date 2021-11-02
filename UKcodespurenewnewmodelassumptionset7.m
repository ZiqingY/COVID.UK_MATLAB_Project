% assumption set 7, almost perfectly match the data, but estimate more deaths on 2nd peak, indicating fatigue almost constant during lockdown easing
set(groot,'defaultLineLineWidth',2.0);
zeta = 1/30;
%zeta = 10;
gamma = 0.4;
sigma = 0.425;
eta = 0.025;
nu = 0.2;
betabar = 3.6*gamma;
betabarv = 5*gamma;
population = 67900000;

seasonalsize = 1;
seasonalposition = 15;

kappa = 1.400e+05;
fatiguesize = 0.575;
fatiguemean = 305;
fatiguesig = 8;

kappa1 = 1.300e+05;
fatiguesize1 = 1.08;
fatiguemean1 = 65;
fatiguesig1 = 25;

kappa2 = 1.05000e+05;
fatiguesize2 = 0.8;
fatiguemean2 = 200;
fatiguesig2 = 200 ;

kappa3 = 1.400e+05;
fatiguesize3 = 0.415;
fatiguemean3 = 300;
fatiguesig3 = 25;

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
tfinal = 570;

y0 = zeros(8,1);
params = zeros(16,1);
y0(1) = S0;
y0(2) = E0;
y0(3) = Ev0;
y0(4) = I0;
y0(5) = Iv0;
y0(6) = R0;
y0(7) = H0;
y0(8) = D0;
params(1) = betabar;
params(2) = kappa;
params(3) = kappa1;
params(4) = kappa2;
params(5) = kappa3;
params(6) = zeta;
params(7) = cases;
params(8) = gamma;
params(9) = sigma;
params(10) = eta;
params(11) = nu;
params(12) = seasonalsize;
params(13) = fatiguesize1;
params(14) = fatiguemean1;
params(15) = fatiguesig1;
params(16) = fatiguesize2;
params(17) = fatiguemean2;
params(18) = fatiguesig2;
params(19) = fatiguesize3;
params(20) = fatiguemean3;
params(21) = fatiguesig3;
params(22) = betabarv;
params(23) = tv;
params(24) = Evbar;
params(25) = seasonalposition;

opts = odeset('Reltol',1e-13,'AbsTol',1e-14,'Stats','on');
[t1,y1] = ode113(@(t,y) bsirSEIRHDassumptionset2(t,y,params),[t0 tfinal],y0,opts);


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

%subplot(2,2,1)
%plot(t1date,D1*population,tdatadate,ukdatanewnew(:,3))
%line([37 37],[0 250000],'linestyle', '--' , 'Color','r', 'LineWidth', 1)
%line([129 129],[0 250000],'linestyle', '--' , 'Color','g', 'LineWidth', 1)
%line([259 259],[0 250000],'linestyle', '--' , 'Color','r', 'LineWidth', 1)
%line([291 291],[0 250000],'linestyle', '--' , 'Color','g', 'LineWidth', 1)
%line([324 324],[0 250000],'linestyle', '--' , 'Color','r', 'LineWidth', 1)
%line([422 422],[0 250000],'linestyle', '--' , 'Color','g', 'LineWidth', 1)
figure(1)
plot(t1date,nu*zeta*H1*population,tdatadate,ukdatanewnew(:,4))
line([37 37],[0 1500],'linestyle', '--' , 'Color','r', 'LineWidth', 1)
line([129 129],[0 1500],'linestyle', '--' , 'Color','g', 'LineWidth', 1)
line([259 259],[0 1500],'linestyle', '--' , 'Color','r', 'LineWidth', 1)
line([291 291],[0 1500],'linestyle', '--' , 'Color','g', 'LineWidth', 1)
line([324 324],[0 1500],'linestyle', '--' , 'Color','r', 'LineWidth', 1)
line([422 422],[0 1500],'linestyle', '--' , 'Color','g', 'LineWidth', 1)

% plot fatigue
t=transpose((1:365))
tnew=transpose((1:570))

kt_1 = (1 + (fatiguesize1-1)*normcdf(t,fatiguemean1,fatiguesig1))*kappa1
kt_2 = (1 + (fatiguesize2-1)*normcdf(t,fatiguemean2,fatiguesig2))*kappa2
kt_3 = (1 + (fatiguesize3-1)*normcdf(t,fatiguemean3,fatiguesig3))*kappa3
k_t = (1 + (fatiguesize-1)*normcdf(t,fatiguemean,fatiguesig))*kappa

%subplot(3,1,2)
%plot(t,k_t, t,kt_1, t,kt_2, t,kt_3)
%line([37 37],[70000 150000],'linestyle', '--' , 'Color','r', 'LineWidth', 1)
%line([129 129],[70000 150000],'linestyle', '--' , 'Color','g', 'LineWidth', 1)
%line([259 259],[70000 150000],'linestyle', '--' , 'Color','r', 'LineWidth', 1)
%line([291 291],[70000 150000],'linestyle', '--' , 'Color','g', 'LineWidth', 1)
%line([324 324],[70000 150000],'linestyle', '--' , 'Color','r', 'LineWidth', 1)
%line([422 422],[70000 150000],'linestyle', '--' , 'Color','g', 'LineWidth', 1)

figure(2)
ktnew = (1 + (fatiguesize1-1)*normcdf(t,fatiguemean1,fatiguesig1))*kappa1.*(t<=129) + (1 + (fatiguesize2-1)*normcdf(t,fatiguemean2,fatiguesig2))*kappa2.*(t>129 & t<=259) + (1 + (fatiguesize3-1)*normcdf(t,fatiguemean3,fatiguesig3))*kappa3.*(t>259)
plot(t,ktnew)
line([37 37],[50000 150000],'linestyle', '--' , 'Color','r', 'LineWidth', 1)
line([129 129],[50000 150000],'linestyle', '--' , 'Color','g', 'LineWidth', 1)
line([259 259],[50000 150000],'linestyle', '--' , 'Color','r', 'LineWidth', 1)
line([291 291],[50000 150000],'linestyle', '--' , 'Color','g', 'LineWidth', 1)
line([324 324],[50000 150000],'linestyle', '--' , 'Color','r', 'LineWidth', 1)
line([422 422],[50000 150000],'linestyle', '--' , 'Color','g', 'LineWidth', 1)