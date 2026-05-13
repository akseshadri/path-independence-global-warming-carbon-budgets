% Fig. 3: Evolution of theta (ratio of warming/CO2-increase)
% (a) Emissions & concentration (b) d(theta)/dt contributions (c) theta vs concentration (d) theta(t)
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

%% taud and taus, tauf
tauD = cd*(beta+nu*gamma)/nu/gamma^2;
taus = cd*(beta+nu*gamma)/beta/gamma;
tauf = cs/(beta+nu*gamma);

%% ODE options
optionsc = odeset('RelTol',1e-2,'AbsTol',1e5,'MaxStep',1.0);
optionst = odeset('RelTol',1e-3,'AbsTol',[1e-4 1e-4],'MaxStep',0.2);

%% Calculate CO2 concentration
[tco2,x2,m1] = co2conc_co2only(gr,tstab,taum2);
mintime = min(tco2); maxtime = 3000; 

%% Calculate emissions
m = diff(m1)./diff(tco2);
t1 = getnpointavg(tco2,2);

%% Solve 2-box model
[s,Y] = ode45(@twoboxmodelco2,[mintime maxtime],[Ts0 Td0],optionst);

Ts = Y(:,1); [Tsmax,indexmax] = max(Ts);

Td = Y(:,2);

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

%% Calculate l1
l1 = 1./taug.*g0./(gtg2-g0);

%% Calculate l2
thetatg2 = interp1(tg,theta,tg2);
l2 = (1./taug).*gtg2./(gtg2-g0).*thetatg2+1/taus*thetatg2;

%% Calculate l3
l3 = (1/tauD+1/taus)*g0./(gtg2-g0).*log(gtg2/g0);

%% Plot Figure 2
xlims = [2000 3000];
maxm = max(m);

i1 = t1 >= min(xlims) & t1 <= max(xlims);
i2 = tg >= min(xlims) & tg <= max(xlims);

figure(2), 
subplot(2,2,1), [hAx,hLine1,hLine2] = plotyy(t1(i1),cfp*m(i1),tg(i2),cfp*g(i2));
xlabel('year'), ylabel(hAx(1),'annual emissions (PgC.year^{-1})') % left y-axis
ylabel(hAx(2),'concentration (PgC)') % right y-axis
grid on

subplot(2,2,4), plot(tg2,thetatg2, 'k','LineWidth',1), xlabel('year'), ylabel('\theta')
xlim(xlims), ylim([0.7 1.4]), hold on, grid on

subplot(2,2,3), plot(cfp*gtg2,thetatg2,'k','LineWidth',1), xlabel('concentration (PgC)'), ylabel('\theta')
hold on, grid on, xlim([645 1300]), %ylim([0.7 1.4])

subplot(2,2,2), plot(tg2,l1-l2,'r','LineWidth',1), hold on, 
plot(tg2,l3,'b','LineWidth',1),  plot(tg2,l1-l2+l3,'k','LineWidth',1),
xlabel('year'), ylabel('contributions to d\theta/dt  (year^{-1})')
xlim(xlims), ylim([-0.006 0.005]), 
legend('l_1 - l_2', 'l_3','d\theta/dt','Location','NorthEast'), legend('boxoff')
grid on

%% Calculate thetaeq
thetaeq = g0./gtg2.*(1+taug/tauD.*log(gtg2/g0));
figure(3), plot(tg2,thetaeq)
    