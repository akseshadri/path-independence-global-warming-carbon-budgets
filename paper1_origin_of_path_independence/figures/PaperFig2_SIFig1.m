% Fig. 2 & SI Fig. 1: Airborne fraction dynamics
% (a) Emissions (b) m1*dr/dt (c) Contributions h1,h2,h3 from Eq.(3-4) (d) Airborne fraction r(t)
%
% Paper: Seshadri (2017), Clim. Dyn. 49:3383-3401, doi:10.1007/s00382-016-3519-3
%

clear, close all

global RFflag gr thist Fhist tstab mu1 mu2 mu3 mu4 tau1 tau2 tau3

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


%% ODE options
optionsc = odeset('RelTol',1e-2,'AbsTol',1e5,'MaxStep',1.0);
optionst = odeset('RelTol',1e-2,'AbsTol',[1e-2 1e-2],'MaxStep',1.0);

%% Load historical RF
[thist,Fhist] = textread('histF.txt');
mintime = min(thist); maxtime = max(thist);

%% Calculate CO2 concentration
[t,g,m1] = co2conc_co2only(gr,tstab,taum2);


%% Calculations for Fig 1

% emissions
m = diff(m1)./diff(t);
t1 = getnpointavg(t,2);

% excess concentrations
eg = g - g(1);

% ratio
r = eg./m1; 

% r'
rp = diff(r)./diff(t);

% m'
mp = diff(m)./diff(t1);
t2 = getnpointavg(t1,2);

% h2
rt1 = interp1(t,r,t1);
h2 = (mu3+mu4-rt1).*m;

% h3
t1z = t1 - min(t1);
dz = t1z(2) - t1z(1);
h3 = NaN(1,1000);
for i = 1:1000
    zi = t1z(1:i);
    t1i = t1z(i);
    mzi = m(1:i);
    h3i = mu3/tau3*exp(-t1i/tau3)*sum(exp(zi/tau3).*mzi.*dz);
    h3(i) = h3i;
end

% f = m1*r'
m1t1 = interp1(t,m1,t1);
f = m1t1.*rp;

% h1
h1 = f(1:1000) - h2(1:1000) + h3;

% h1a, approx to h1
h1a = NaN(1,1000);
for i = 1:1000
    zi = t1z(1:i);
    t1i = t1z(i);
    mzi = m(1:i);
    mti = m(i);
    h1ai = mu2*(mti - 1/tau2*exp(-t1i/tau2)*sum(exp(zi/tau2).*mzi.*dz));% + mu1*(mti - 1/tau1*exp(-t1i/tau1)*sum(exp(zi/tau1).*mzi.*dz));
    h1a(i) = h1ai;
end

%% Limits
xlims = [1850 2400];
i = t >= min(xlims) & t <= max(xlims);
maxm = max(m);
maxr = max(r(i)); 
minr = min(r(i));

%% Graphs of Fig 1
figure(1), 
subplot(2,2,1), plot(t1,m*cfp,'k','LineWidth',1), xlabel('year'), ylabel('emissions (PgC.year^{-1})')
xlim(xlims), %ylim([0 1.02*maxm]), 
hold on, grid on

subplot(2,2,4), plot(t,r, 'k','LineWidth',1), xlabel('year'), ylabel('airborne fraction')
xlim(xlims), ylim([minr maxr]), hold on, grid on

subplot(2,2,2), plot(t1,f*cfp,'k','LineWidth',1), xlabel('year'), ylabel('m_1*dr/dt (PgC.year^{-1})')
xlim(xlims), hold on, %ylim([1.02*min(f) 1.02*max(f)]), 
grid on

subplot(2,2,3), plot(t1(1:1000),h1*cfp,'r','LineWidth',1), hold on, plot(t1,h2*cfp,'g','LineWidth',1), 
plot(t1(1:1000),-h3*cfp,'b','LineWidth',1), 
xlabel('year'), ylabel('contributions to m_1*dr/dt (PgC.year^{-1})')
xlim(xlims), legend('  h_1', '  h_2','- h_3','Location','NorthEast'), legend('boxoff')
%ylim([1.02*min(h2) 1.02*max(h1)]), 
grid on
plot(t1(1:1000),h1a*cfp,'r:','LineWidth',1)

%% Calculations for Fig 2
m2 = cumsum(m1); 
m3 = cumsum(m2);
m4 = cumsum(m3);
m5 = cumsum(m4);
m6 = cumsum(m5);

l1 = mu3/tau3*(m1-m2/tau3);
l2 = mu3/tau3^3*(m3-m4/tau3);
l3 = mu3/tau3^5*(m5-m6/tau3);

%% Graphs of Fig 2
  
figure(2), 
subplot(2,1,1), plot(t1,m*cfp,'k','LineWidth',1.5), xlabel('year'), ylabel('emissions (PgC.year^{-1})'), xlim(xlims), %ylim([0 1.02*max(m)])
subplot(2,1,2), plot(t,l1*cfp,'r','LineWidth',1.5), hold on, plot(t,l2*cfp,'g','LineWidth',1.5), plot(t,l3*cfp,'b','LineWidth',1.5), plot(t1(1:1000),h3*cfp,'k','LineWidth',1.5),
xlabel('year'), ylabel('first few terms in series approximating h_3 (PgC.year^{-1})'), 
legend('h_3_1','h_3_2','h_3_3','h_3','Location','NorthWest'), legend('boxoff'), xlim(xlims), %ylim([-0.01 1.02*max(h3)])

%% Output for discussion
[mindyr,itaur] = min(abs(t1-2015));
ri = r(itaur);
m1i = m1(itaur)*cfp;
mi = m(itaur)*cfp;

taur1 = -m1i*ri/cfp/h1(itaur);
taur2 = -m1i*ri/cfp/h2(itaur);
taur3 = m1i*ri/cfp/h3(itaur);

taur = 1/(1/taur1 + 1/taur2 + 1/taur3);