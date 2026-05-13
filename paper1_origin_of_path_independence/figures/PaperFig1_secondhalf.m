% Fig. 1 (continued): Additional emissions scenarios for EBM verification
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
    
    %% Calculate parameters of EBM
    b1 = (beta + nu*gamma) / cs;
    b2 = nu*gamma / cs;
    b3 = gamma / cs;
    b4 = gamma / cs;
    
    %% taud and tau1, tau2
    tauD = cd*(beta+nu*gamma)/nu/gamma^2;
    tau1 = cd*(beta+nu*gamma)/beta/gamma;
    tau2 = cs/(beta+nu*gamma);
    
    %% Simulate approx
    slist = s; 
    slist_zeroed = slist - min(slist);
    
    sl = linspace(min(slist),max(slist),5000);
    Fl = interp1(tco2,F,sl);       
        
    Tsest = NaN(numel(sl)-1,1);
    Tdest = NaN(numel(sl)-1,1);
    
    for i  = 1:numel(sl)-1
        si = sl(i);
        zi = sl(1:i); dzi = diff(sl(1:i+1));
        Fllisti = Fl(1:i);
        Fli = Fl(i);
        
        Tsesti = tau2/cs*(Fli + 1/tauD*exp(-si/tau1)*sum(exp(zi/tau1).*Fllisti.*dzi));
        Tdesti = epsilon*b3/cs/b1*(exp(-si/tau1)*sum(exp(zi/tau1).*Fllisti.*dzi));
        
        Tsest(i) = Tsesti;
        Tdest(i) = Tdesti;
    end
    
       
    %% Plot
    figure(1),
    subplot(2,2,1),
    plot(t1,m*cfp, char(pal(maincount)), 'LineWidth', 1.0), hold on, ...
        xlabel('year'), ylabel('emissions (PgC.year^{-1})'),...
        xlim([1750 2500])
    subplot(2,2,2), 
    plot(tco2,m1*cfp, char(pal(maincount)), 'LineWidth', 1.0), hold on, ...
        xlabel('year'), ylabel('cumulative emissions (PgC)'),...
        xlim([1750 2500])
    subplot(2,2,3),
    plot(t1,Ft1,char(pal(maincount)), 'LineWidth', 1.0), hold on, ...
        xlabel('year'), ylabel('radiative forcing, F(t)  (W m^{-2})'),...
        xlim([1750 2500])
    
    subplot(2,2,4)
    plot(s,Ts, char(palc(maincount)), 'LineWidth', 1.0), hold on, ...
        plot(sl(1:numel(sl)-1),Tsest, char(pala(maincount)), 'LineWidth', 1.0), xlabel('year'), ylabel('global warming, T_s  (K)'),...
        xlim([1750 2500])
    
    
end





