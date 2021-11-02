%% Codes in this file are for evaluation of UK lockdowns

t = transpose((1:365))
tnew=transpose((1:570))

seasonalsize = 1;
seasonalposition = 15
psi = seasonalsize*(cos((tnew+seasonalposition)*2*pi/365)-1)/2;

kappa1 = 1.300e+05;
fatiguesize1 = 1.08;
fatiguemean1 = 70;
fatiguesig1 = 36;

kappa2 = 1.05000e+05;
fatiguesize2 = 0.8;
fatiguemean2 = 200;
fatiguesig2 = 100 ;

kappa3 = 1.500e+05;
fatiguesize3 = 0.415;
fatiguemean3 = 290;
fatiguesig3 = 25;

% construct the fatigue parameter
knew = (1 + (fatiguesize1-1)*normcdf(tnew,fatiguemean1,fatiguesig1))*kappa1.*(tnew<=129) + (1 + (fatiguesize2-1)*normcdf(tnew,fatiguemean2,fatiguesig2))*kappa2.*(tnew>129 & tnew<=259) + (1 + (fatiguesize3-1)*normcdf(tnew,fatiguemean3,fatiguesig3))*kappa3.*(tnew>259)
format long
load t1H1forHtasp7.mat
Ht = t1H1forHtasp7((1:570))
kneweffect = -knew.*Ht*0.2/30
t1datenew = datetime('15-Feb-2020')+tnew
tdatadatenew = datetime('15-Feb-2020')+tnew

% convert the dimension of 'hospitalization' 
figure(1)
subplot(2,1,1)
plot(t1datenew, Ht)
subplot(2,1,2)
plot(t1date, H1)


% compare seasonality and fatigue parameters
figure(2)
subplot(2,2,1)
plot(t1datenew, psi)
title('seasonality','FontSize',16)
line([37 37],[-1.2 0.2],'linestyle', '--' , 'Color','r', 'LineWidth', 1)
line([129 129],[-1.2 0.2],'linestyle', '--' , 'Color','g', 'LineWidth', 1)
line([259 259],[-1.2 0.2],'linestyle', '--' , 'Color','r', 'LineWidth', 1)
line([291 291],[-1.2 0.2],'linestyle', '--' , 'Color','g', 'LineWidth', 1)
line([324 324],[-1.2 0.2],'linestyle', '--' , 'Color','r', 'LineWidth', 1)
line([422 422],[-1.2 0.2],'linestyle', '--' , 'Color','g', 'LineWidth', 1)


subplot(2,2,2)
plot(t1datenew, knew/100000-1.5, t1datenew, psi)
title('Compare fatigue parameter k and seasonality', 'FontSize',16)
line([37 37],[-2 0.5],'linestyle', '--' , 'Color','r', 'LineWidth', 1)
line([129 129],[-2 0.5],'linestyle', '--' , 'Color','g', 'LineWidth', 1)
line([259 259],[-2 0.5],'linestyle', '--' , 'Color','r', 'LineWidth', 1)
line([291 291],[-2 0.5],'linestyle', '--' , 'Color','g', 'LineWidth', 1)
line([324 324],[-2 0.5],'linestyle', '--' , 'Color','r', 'LineWidth', 1)
line([422 422],[-2 0.5],'linestyle', '--' , 'Color','g', 'LineWidth', 1)


subplot(2,2,3)
plot(t1datenew, psi, t1datenew, kneweffect, t1datenew, psi + kneweffect)
title('Compare effects of seasonality and fatigue parameter k', 'FontSize',16)
line([37 37],[-2.5 0.5],'linestyle', '--' , 'Color','r', 'LineWidth', 1)
line([129 129],[-2.5 0.5],'linestyle', '--' , 'Color','g', 'LineWidth', 1)
line([259 259],[-2.5 0.5],'linestyle', '--' , 'Color','r', 'LineWidth', 1)
line([291 291],[-2.5 0.5],'linestyle', '--' , 'Color','g', 'LineWidth', 1)
line([324 324],[-2.5 0.5],'linestyle', '--' , 'Color','r', 'LineWidth', 1)
line([422 422],[-2.5 0.5],'linestyle', '--' , 'Color','g', 'LineWidth', 1)
legend('psi', 'keffect', 'psi+keffect')
legend('boxoff')


subplot(2,2,4)
beta570 = betabar*exp(-kappa*nu*zeta*Ht+psi)
betav570 = betabarv*exp(-kappa*nu*zeta*Ht+psi)
plot(t1datenew, beta570, t1datenew, betav570, t1datenew, kneweffect, t1datenew, psi)
title('beta, effect of psi, effect of k', 'FontSize',16)
line([37 37],[-2 1.5],'linestyle', '--' , 'Color','r', 'LineWidth', 1)
line([129 129],[-2 1.5],'linestyle', '--' , 'Color','g', 'LineWidth', 1)
line([259 259],[-2 1.5],'linestyle', '--' , 'Color','r', 'LineWidth', 1)
line([291 291],[-2 1.5],'linestyle', '--' , 'Color','g', 'LineWidth', 1)
line([324 324],[-2 1.5],'linestyle', '--' , 'Color','r', 'LineWidth', 1)
line([422 422],[-2 1.5],'linestyle', '--' , 'Color','g', 'LineWidth', 1)
legend('beta', 'k', 'psi')
legend('boxoff')


% plot the simulations for accumulated death and daily death
figure(3)
subplot(2,1,1)
plot(t1date,D1*population,tdatadate,ukdatanewnew(:,3))
line([37 37],[0 250000],'linestyle', '--' , 'Color','r', 'LineWidth', 1)
line([129 129],[0 250000],'linestyle', '--' , 'Color','g', 'LineWidth', 1)
line([259 259],[0 250000],'linestyle', '--' , 'Color','r', 'LineWidth', 1)
line([291 291],[0 250000],'linestyle', '--' , 'Color','g', 'LineWidth', 1)
line([324 324],[0 250000],'linestyle', '--' , 'Color','r', 'LineWidth', 1)
line([422 422],[0 250000],'linestyle', '--' , 'Color','g', 'LineWidth', 1)
subplot(2,1,2)
plot(t1date,nu*zeta*H1*population,tdatadate,ukdatanewnew(:,4))
line([37 37],[0 1500],'linestyle', '--' , 'Color','r', 'LineWidth', 1)
line([129 129],[0 1500],'linestyle', '--' , 'Color','g', 'LineWidth', 1)
line([259 259],[0 1500],'linestyle', '--' , 'Color','r', 'LineWidth', 1)
line([291 291],[0 1500],'linestyle', '--' , 'Color','g', 'LineWidth', 1)
line([324 324],[0 1500],'linestyle', '--' , 'Color','r', 'LineWidth', 1)
line([422 422],[0 1500],'linestyle', '--' , 'Color','g', 'LineWidth', 1)
