% Figs. 5-6 variant: Second set of emissions scenarios
%
% Paper: Seshadri (2017), Clim. Dyn. 49:3383-3401, doi:10.1007/s00382-016-3519-3
%


clear

global RFflag gr cs cd beta nu gamma tco2 x2 tstab mu1 mu2 mu3 mu4 tau1 tau2 tau3 x2_PI tn1 tn2 sn

RFflag = 2;

pal = {'k','b-.','r--'};
palc = {'k','b','r'};
pala = {'k--','b--','r--'};

taum2list = [40 42 40];

%% conversion factor
cfp = 1/(3.67*1e9*120/9.5e14)/1e3; 

%% economy parameters
grlist = [0.03 0.03 0.03];
tstablist = [60 120 60];

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

%% noise inputs
tn1list = [120 120 120];
tn2list = [450 450 450];
snlist = [0 0 1];

for maincount = 1:3
    
     %% emission parameters
    gr = grlist(maincount);
    tstab = tstablist(maincount);
    taum2 = taum2list(maincount);
    
    %% Noise parameters
    tn1 = tn1list(maincount);
    tn2 = tn2list(maincount);
    sn = snlist(maincount);
    
    %% Calculate CO2 concentration
    [tco2,x2,m1] = co2conc_co2onlywp2(gr,tstab,taum2);
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
    
    %% interpolate to s
    m1s = interp1(tco2,m1,s);
    
    %% Calculate errors and contributions to dx
    dz = t1(2) - t1(1);
    
    err_r = NaN(numel(t1),1);
    err_t = NaN(numel(t1),1);
    Tsapprox = NaN(numel(t1),1);
    Tsapproxf = NaN(numel(t1),1);
    Tsapproxs = NaN(numel(t1),1);
    
    u_m1 = NaN(numel(t1),1);
    u_r = NaN(numel(t1),1);
    u_t = NaN(numel(t1),1);
    
    
    for i = 1:numel(t1)
        taum1i = taum1(i);
        tauri = taur(i);
        gbym1i = 1./m1byg(i);
        gbyri = 1./rbyg(i);
        gi = gt1(i);
        m1i = m1t1(i);
        mi = m(i);
        
        ri = rt1(i);
        
        intF = sum(Ft1(1:i));
        Fi = Ft1(i);
        
        Tsapproxi = Tst1(i);
        Tsapproxfi = tauf/cs*(Fi);
        Tsapproxsi = Tsapproxi - Tsapproxfi;
        
        intm1bygi = sum(m1byg(1:i));
        intrbygi = sum(rbyg(1:i));
        
        er_ri = taum1i/tauri*(1+1/tauD*gi/m1i*intm1bygi*dz)/(1+1/tauD*gi/ri*intrbygi*dz);
        er_ti = taum1i/tauD/(1+(gi-g0)/gi/log(gi/g0)*(1+1/tauD*gi/ri*intrbygi*dz));
        
        err_r(i) = abs(er_ri);
        err_t(i) = abs(er_ti);
        Tsapprox(i) = Tsapproxi;
        Tsapproxf(i) = Tsapproxfi;
        Tsapproxs(i) = Tsapproxsi;
        
        u_m1(i) = (ri/gi + 1/tauD*intrbygi*dz);
        u_r(i) = (m1i/gi + 1/tauD*intm1bygi*dz)/mi*dr(i);
        u_t(i) = (1/tauD*log(gi/g0))/mi; 
    end
    
    x_m1 = 5.35*tauf/cs*u_m1;
    x_r = 5.35*tauf/cs*u_r;
    x_t = 5.35*tauf/cs*u_t;
    
    x_tot = x_m1+x_r+x_t;
    
    %% Calculate dT/dm1
    dTst1 = diff(Tst1)./diff(m1t1);
    m1t12 = getnpointavg(m1t1,2);
    
    %% Plot Figure 6
    
   xlims = [0 1400*cfp];
    
   figure(7),
    subplot(3,2,1), plot(m1t1*cfp,m*cfp,char(pal(maincount)), 'LineWidth',1.5), hold on, ...
        xlabel('cumulative CO_2 emissions (PgC)'), ylabel('emissions (PgC.year^{-1})'), ...
        xlim(xlims), ylim([0 20]),...
        grid on,
    
    subplot(3,2,2), plot(m1t12*cfp,dTst1/cfp,char(pal(maincount)), 'LineWidth',1.5), hold on,...
        xlabel('cumulative CO_2 emissions (PgC)'), ylabel('sensitivitiy of global warming dx/dm_1 (K.PgC)'), grid on,...
        ylim([0 3e-3]), xlim(xlims), 
    
    subplot(3,2,4), semilogy(m1t1*cfp,abs(taum1./tauD),char(pal(maincount)),'LineWidth',1.5),grid on, hold on,...
        xlabel('cumulative CO_2 emissions (PgC)'), ylabel('abs(\varsigma_m_1 / \tau_D)'), ylim([0.01 10]), xlim(xlims), 
    
    
    subplot(3,2,3), semilogy(m1t1*cfp,abs(taum1./taur),char(pal(maincount)),'LineWidth',1.5),grid on, hold on,...
        xlabel('cumulative CO_2 emissions (PgC)'), ylabel('abs(\varsigma_r / \varsigma_m_1)'), ylim([0.01 10]), xlim(xlims), 
    
    
    subplot(3,2,5), plot(m1t1*cfp,rt1,char(pal(maincount)), 'LineWidth',1.5), hold on,...
        xlabel('cumulative CO_2 emissions (PgC)'), ylabel('airborne fraction r(t)'), grid on, xlim(xlims), 
    ylim([0.3 0.6])
    
    subplot(3,2,6), plot(m1*cfp,theta,char(pal(maincount)), 'LineWidth',1.5), hold on,...
        xlabel('cumulative CO_2 emissions (PgC)'), ylabel('global warming/carbon increase \theta(t)'), grid on, xlim(xlims), 
    ylim([0.6 1.1])
    
    figure(8),
    plot(m1t1*cfp,Tst1,char(pal(maincount)), 'LineWidth',1.5), hold on, ...
        xlabel('cumulative CO_2 emissions (PgC)'), ylabel('global warming (K)'), ...
        xlim([0 4000]), ylim([0 5])
    
    figure(9),
    plot(t1,rt1,char(pal(maincount)), 'LineWidth',1.5), hold on, ...
        xlabel('year'), ylabel('airborne fraction')
    
end





