% Fig. 4: Timescales for path independence condition (Eq. 13)
% Compares zeta_r(t), zeta_m1(t), and tau_D. Requires: zeta_r >> zeta_m1 << tau_D
%
% Paper: Seshadri (2017), Clim. Dyn. 49:3383-3401, doi:10.1007/s00382-016-3519-3
%

clear, close all

global RFflag gr cs cd beta nu gamma tco2 x2 tstab mu1 mu2 mu3 mu4 tau1 tau2 tau3 x2_PI

RFflag = 2;

palexact = {'k','b','r'};
palapprox = {'k:','b:','r:'};

taum2 = 40;

%% conversion factor
cfp = 1/(3.67*1e9*120/9.5e14)/1e3;

%% economy parameters
gr = 0.03;
tstab = 60;

%% carbon cycle model parameters
mu1 = 0.276;
mu2 = 0.282;
mu3 = 0.224;

mu4 = 1 - (mu1+mu2+mu3);
tau1 = 4.3;
tau2 = 36.5;
tau3 = 394;
x2_PI = 285; % ppm

%% 2-box model parameters
cs0 = 8.2; % W*yr / m^2.K
cd0 = 109;
ctot0 = cs0 + cd0;
cs = cs0; cd = cd0;
epsilon = 0.075; r = 1/epsilon;

gamma = 0.67; % W/m^2.K
nu = 1.28; % efficacy of ocean heat uptake
DeltaTECS = 3.1; % K

beta = 5.35*log(2)/DeltaTECS; % W/m^2.K

Ts0 = 0.0; Td0 = 0.0;

%% ODE options
optionsc = odeset('RelTol',1e-2,'AbsTol',1e5,'MaxStep',1.0);
optionst = odeset('RelTol',1e-3,'AbsTol',[1e-4 1e-4],'MaxStep',0.2);

%% Calculate CO2 concentration
[tco2,x2,m1] = co2conc_co2only(gr,tstab,taum2);
mintime = min(tco2); maxtime = 3000; 

%% Calculate emissions
m = diff(m1)./diff(tco2);
t1 = getnpointavg(tco2,2);

%% Calculate forcing
F = 5.35*log(x2/x2_PI);
Ft1 = interp1(tco2,F,t1);

%% Solve 2-box model
[s,Y] = ode45(@twoboxmodelco2,[mintime maxtime],[Ts0 Td0],optionst);

Ts = Y(:,1); [Tsmax,indexmax] = max(Ts);

Td = Y(:,2);

%% Calculate u
Hx = gamma;
Hy = -gamma;
Bx = beta + (nu - 1)*gamma;
By = - (nu - 1)*gamma;

%% taud and taus, tauf
tauD = cd*(beta+nu*gamma)/nu/gamma^2;
taus = cd*(beta+nu*gamma)/beta/gamma;
tauf = cs/(beta+nu*gamma);

%% Calculate u and theta
g0 = x2_PI;

tg = min(s):3000;
g = interp1(tco2,x2,tg);

u = cs/5.35/tauf*interp1(s,Ts,tg);

theta = u.*g0./(g-g0);

%% Calculate g' and taug
gp = diff(g)./diff(tg);
tg2 = getnpointavg(tg,2);

gtg2 = interp1(tg,g,tg2);

taug = gtg2./gp;

%% Calculate theta'
thetap = diff(theta)./diff(tg);

%% Calculate r
r = (g-g0)./m1; 

%% Calculate taum1 and taur
m1t1 = interp1(tco2,m1,t1);
taum1 = m1t1./m;

dr = diff(r)./diff(tco2);
rt1 = interp1(tco2,r,t1);
taur = -rt1./dr;

%% Interpolate to t1
gt1 = interp1(tco2,g,t1);

%% Calculate m1/g and r/g
m1byg = m1t1./gt1;
rbyg = rt1./gt1;

%% Calculate global warming
Tst1 = interp1(s,Ts,t1);

%% Calculate errors
dz = t1(2) - t1(1);

err_r = NaN(numel(t1),1);
err_t = NaN(numel(t1),1);
Tsapprox = NaN(numel(t1),1);

for i = 1:numel(t1)
    taum1i = taum1(i);
    tauri = taur(i);
    gbym1i = 1./m1byg(i);
    gbyri = 1./rbyg(i);
    gi = gt1(i);
    m1i = m1t1(i);
    ri = rt1(i);
    
    intF = sum(Ft1(1:i));
    Fi = Ft1(i);
    
    Tsapproxi = tauf/cs*(Fi+1/tauD*intF*dz);
    
    intm1bygi = sum(m1byg(1:i));
    intrbygi = sum(rbyg(1:i));
    
    er_ri = taum1i/tauri*(1+1/tauD*gi/m1i*intm1bygi*dz)/(1+1/tauD*gi/ri*intrbygi*dz);
    er_ti = taum1i/tauD/(1+(gi-g0)/gi/log(gi/g0)*(1+1/tauD*gi/ri*intrbygi*dz));
    
    err_r(i) = abs(er_ri);
    err_t(i) = abs(er_ti);
    Tsapprox(i) = Tsapproxi;
end

%% Plot Figures 3
% fig 3
figure(3), 
subplot(2,1,1), plot(t1,m*cfp, 'k'), xlabel('year'), ylabel('annual emissions (PgC.year^{-1})'), xlim([1750 2200]), grid on

subplot(2,1,2), semilogy(t1,abs(taur),'r', 'LineWidth',1.0), hold on, semilogy(t1,ones(size(t1))*tauD,'k', 'LineWidth',1.0), semilogy(t1,taum1,'b', 'LineWidth',1.0),...
    xlabel('year'), ylabel('timescale (years)'), legend('abs(\varsigma_r)','\tau_D','\varsigma_m_1'),...
    xlim([1750 2200]), ylim([0 10000]), legend boxoff, grid on

n =10;
tauravg = getnpointavg(abs(taur),10);
t15 = getnpointavg(t1,10);
taum1avg = getnpointavg(taum1,10);

plot(t15,tauravg,'r-.','LineWidth',1.0), plot(t15,taum1avg,'b-.','LineWidth',1.0) 

%% Figure for presentation
figure(1), 

plot(t1,abs(taur),'r', 'LineWidth',2.0), hold on, plot(t1,taum1,'b', 'LineWidth',2.0),...
    xlabel('year'), ylabel('timescale (years)'), legend('airborne fraction','cumulative CO_2 emissions'),...
    xlim([2000 2100]), ylim([10 1000]), grid on, 




    