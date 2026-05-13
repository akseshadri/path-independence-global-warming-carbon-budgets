% SI Fig. 4: Airborne fraction differences and their contributions
%
% Paper: Seshadri (2017), Clim. Dyn. 49:3383-3401, doi:10.1007/s00382-016-3519-3
%

clear, close all

global RFflag gr thist Fhist tstab mu1 mu2 mu3 mu4 tau1 tau2 tau3 tn1 tn2 sn

RFflag = 2;

pal = {'k','b-.','r--'};
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


%% ODE options
optionsc = odeset('RelTol',1e-2,'AbsTol',1e5,'MaxStep',1.0);
optionst = odeset('RelTol',1e-2,'AbsTol',[1e-2 1e-2],'MaxStep',1.0);

%% noise inputs
tn1list = [120 120 120];
tn2list = [450 450 450];
snlist = [0 0 1];

%% Load historical RF
[thist,Fhist] = textread('histF.txt');
mintime = min(thist); maxtime = max(thist);

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
    [t,g,m1] = co2conc_co2onlywp2(gr,tstab,taum2);
    
    
    %% Calculations for Fig 4
    
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
    
    n = 1200;
    
    % h3
    t1z = t1 - min(t1);
    dz = t1z(2) - t1z(1);
    h3 = NaN(1,n);
    for i = 1:n
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
    h1 = f(1:n) - h2(1:n) + h3;
    
    % h1a, approx to h1
    h1a = NaN(1,n);
    for i = 1:n
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
    xlims = [0 1400*cfp]; ylims = [-1.5*cfp 1*cfp];
    
    figure(4),
    subplot(2,2,1), plot(cfp*m1t1(1:n),h1*cfp,char(pal(maincount)),'LineWidth',1), hold on, ...
        xlabel('cumulative CO_2 emissions (PgC)'), ylabel('h_1 (PgC.year^{-1})'),...
        xlim(xlims), ylim(ylims), grid on
    
    subplot(2,2,2), plot(cfp*m1t1,h2*cfp,char(pal(maincount)),'LineWidth',1), hold on, ...
        xlabel('cumulative CO_2 emissions (PgC)'), ylabel('h_2 (PgC.year^{-1})'),...
        xlim(xlims), ylim(ylims), grid on
    
    subplot(2,2,3), plot(cfp*m1t1(1:n),-h3*cfp,char(pal(maincount)),'LineWidth',1), hold on, ...
        xlabel('cumulative CO_2 emissions (PgC)'), ylabel('h_3 (PgC.year^{-1})'),...
        xlim(xlims), ylim(ylims), grid on
    
    subplot(2,2,4), plot(cfp*m1t1,f*cfp,char(pal(maincount)),'LineWidth',1), xlabel('cumulative CO_2 emissions (PgC)'), ylabel('m_1*dr/dt (PgC.year^{-1})'),...
    xlim(xlims), ylim(ylims), hold on,  grid on

       
end
