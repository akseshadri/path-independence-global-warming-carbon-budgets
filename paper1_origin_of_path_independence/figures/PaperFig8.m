% Fig. 8: Effect of f_long (fraction from long carbon cycle time-constants)
% Shows slow carbon cycle time-constants are *essential* for path independence.
% If f_long is too small, the path independent warming-vs-cumulative-emissions relation breaks down.
%
% Paper: Seshadri (2017), Clim. Dyn. 49:3383-3401, doi:10.1007/s00382-016-3519-3
%

clear, close all

global RFflag gr cs cd beta nu gamma tco2 x2 tstab mu1 mu2 mu3 mu4 tau1 tau2 tau3 x2_PI

RFflag = 2;

pal1 = {'g','m','c'};
pal2 = {'g--','m--','c--'};
pal3 = {'g:','m:','c:'};

palf = {'k--','b--','r--'};
pals = {'k-.','b-.','r-.'};  
palds = {'k:','b:','r:'};
palpoint = {'k.','b.','r.'};

taum2list = [60 40 inf];
fracstlist = [0.276+0.282 0.75 1.0];

%% conversion factor
cfp = 1/(3.67*1e9*120/9.5e14)/1e3;

%% economy parameters
grlist = [0.015 0.04 0];
tstablist = [200 50 250];

%% 2-box model parameters
cs0 = 4 * 3.48 / 1.5; % W*yr / m^2.K

for fracstcount = 1:3
    fracst = fracstlist(fracstcount);
    
    mu10 = 0.276; mu20 = 0.282; mu30 = 0.224; mu40 = 1-mu10-mu20-mu30;
    
    %% carbon cycle model parameters
    mu1 = fracst*mu10/(mu10+mu20);
    mu2 = fracst*mu20/(mu10+mu20);
    mu3 = (1-fracst)*mu30/(mu30+mu40);
    mu4 = (1-fracst)*mu40/(mu30+mu40);
    tau1 = 4.3;
    tau2 = 36.5;
    tau3 = 394;
    x2_PI = 285; % ppm
    
    %% Climate model parameters
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
    
    
    for maincount = 1:3
        
        %% emission parameters
        gr = grlist(maincount);
        tstab = tstablist(maincount);
        taum2 = taum2list(maincount);
        
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
        taurat_t = NaN(numel(t1),1);
        sf_t = NaN(numel(t1),1);
               
        for i = 1:numel(t1)
            taum1i = taum1(i);
            tauri = taur(i);
            m1bygi = m1byg(1:i);
            rbygi = rbyg(1:i);
            
            gi = gt1(i);
            m1i = m1t1(i);
            mi = m(i);
            
            ri = rt1(i);
            
            zi = t1(1:i);
            ti = t1(i);
            
            errn = 1+1/tauD*gi/m1i*exp(-ti/taus)*sum(exp(zi/taus).*m1bygi*dz);
            errd = 1+1/tauD*gi/ri*exp(-ti/taus)*sum(exp(zi/taus).*rbygi*dz);
            er_ri = taum1i/tauri*errn/errd;
            
            err_r(i) = abs(er_ri);
            taurat_t(i) = taum1i/tauri;
            sf_t(i) = errn/errd;
        end
        
        
        %% Plot Figure 6
        if fracstcount == 1
            pal = pal1;
        elseif fracstcount == 2
            pal = pal2;
        else
            pal = pal3;
        end
        
        figure(9),
        subplot(2,2,1), semilogy(t1,err_r,char(pal(maincount)), 'LineWidth',1.5), hold on, ...
            xlabel('year'), ylabel('relative error from airborne fraction, abs(\Deltau_r/\Deltau_m_1)'), ...
            xlim([1900 2200]), ylim([1e-2 1])
        
        i = t1 <= 2200;
        grid on,
        subplot(2,2,2),
        plot(cfp*m1t1(i),Tst1(i),char(pal(maincount)), 'LineWidth',1.5), hold on, ...
            xlabel('cumulative CO_2 emissions (PgC)'), ylabel('global warming from CO_2 (K)'), ...
            xlim([0 1000*cfp]), ylim([0 3]), grid on
        
          %% SI fig
        figure(9),
        subplot(2,2,3), semilogy(t1,taurat_t,char(pal(maincount)), 'LineWidth',1.5), hold on, ...
            xlabel('year'), ylabel('ratio of timescales, abs(\varsigma_r / \tau_D)'), ...
            xlim([1900 2200]), ylim([1e-2 10]), grid on
        
        subplot(2,2,4),
         semilogy(t1,sf_t,char(pal(maincount)), 'LineWidth',1.5), hold on, ...
            xlabel('year'), ylabel('second factor in eq. (11), f_r'), ...
            xlim([1900 2200]), ylim([1e-1 10]), grid on
        
    end
end




