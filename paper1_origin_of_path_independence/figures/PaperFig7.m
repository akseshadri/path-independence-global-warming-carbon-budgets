% Fig. 7: Effect of damping timescale tau_D on path independence
% Shows path independence holds even for unrealistically short tau_D.
% Large deep-ocean heat capacity is NOT essential for path independence.
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
ratclist = [1/0.075 5 2];

%% conversion factor
cfp = 1/(3.67*1e9*120/9.5e14)/1e3;

%% economy parameters
grlist = [0.015 0.04 0];
tstablist = [200 50 250];

%% carbon cycle model parameters
mu1 = 0.276*1.0;
mu2 = 0.282*1.0;
mu3 = 0.224*1.0; %0.224
mu4 = 1 - (mu1+mu2+mu3);
tau1 = 4.3;
tau2 = 36.5;
tau3 = 394;
x2_PI = 285; % ppm

%% 2-box model parameters
cs0 = 8.2; % W*yr / m^2.K

for ratccount = 1:3
    ratc = ratclist(ratccount);
    cd0 = ratc*cs0;
    cs = cs0;
    cd = cd0;
          
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
        epsilon = 1/ratc;
       
               
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
        
        %% Calculate r/g and g/g0
        rbyg = rt1./gt1;
        gbyg0 = gt1/g0;
        
        %% Calculate global warming
        Tst1 = interp1(s,Ts,t1);
        
        %% interpolate to s
        m1s = interp1(tco2,m1,s);
        
        %% Calculate errors and contributions to dx
        dz = t1(2) - t1(1);
        
        err_t = NaN(numel(t1),1);
        taurat_t = NaN(numel(t1),1);
        sf_t = NaN(numel(t1),1);
        
        for i = 1:numel(t1)
            taum1i = taum1(i);
            tauri = taur(i);
            
            gi = gt1(i);
            m1i = m1t1(i);
           
            rbygi = rbyg(1:i);
            gbyg0i = gbyg0(1:i);
            
            zi = t1(1:i);
            ti = t1(i);
            dz = t1(2) - t1(1);
            
            ertn = log(gbyg0(i)) - 1/taus*exp(-ti/taus)*sum(exp(zi/taus).*log(gbyg0i)*dz);
            ertd = (gi - g0)/gi + m1i/tauD*exp(-ti/taus)*sum(exp(zi/taus).*rbygi*dz);
            
            er_ti = taum1i/tauD*ertn/ertd;
                     
            
            err_t(i) = abs(er_ti);
            sf_t(i) = ertn/ertd;
            taurat_t(i) = taum1i/tauD;
            
        end
        
        
        %% Plot Figure 
        if ratccount == 1
            pal = pal1;
        elseif ratccount == 2
            pal = pal2;
        else
            pal = pal3;
        end
        
            
        figure(9),
        
        subplot(3,2,1),
        plot(t1,m*cfp, char(pal(maincount)), 'LineWidth', 1.0), hold on, ...
            xlabel('year'), ylabel('emissions (PgC.year^{-1})'),...
            xlim([1900 2200]), grid on
        subplot(3,2,2),
        plot(tco2,m1*cfp, char(pal(maincount)), 'LineWidth', 1.0), hold on, ...
            xlabel('year'), ylabel('cumulative emissions (PgC)'),...
            xlim([1900 2200]), grid on
        
    
        subplot(3,2,3), semilogy(t1,err_t,char(pal(maincount)), 'LineWidth',1.5), hold on, ...
            xlabel('year'), ylabel('relative error from time, abs(\Deltau_t/\Deltau_m_1)'), ...
            xlim([1900 2200]), ylim([1e-2 1])
        
        i = t1 <= 2200;
        grid on,
        subplot(3,2,4),
        plot(cfp*m1t1(i),Tst1(i),char(pal(maincount)), 'LineWidth',1.0), hold on, ...
            xlabel('cumulative CO_2 emissions (PgC)'), ylabel('global warming from CO_2 (K)'), ...
            xlim([0 cfp*1000]), ylim([0 4.0]), grid on
        
        %figure(9),
        %subplot(3,2,1), semilogy(t1,err_t,char(pal(maincount)), 'LineWidth',1.5), hold on, ...
        %    xlabel('year'), ylabel('relative error from time, abs(\Deltau_t/\Deltau_m_1)'), ...
        %    xlim([1900 2200]), ylim([1e-2 1])
        
        
        subplot(3,2,5), semilogy(t1,taurat_t,char(pal(maincount)), 'LineWidth',1.5), hold on, ...
            xlabel('year'), ylabel('ratio of timescales, abs(\varsigma_m_1 / \tau_D)'), ...
            xlim([1900 2200]), ylim([1e-1 10]), grid on
        
        
        subplot(3,2,6),
        semilogy(t1,sf_t,char(pal(maincount)), 'LineWidth',1.5), hold on, ...
            xlabel('year'), ylabel('second factor in eq. (12), f_t'), ...
            xlim([1900 2200]), ylim([1e-1 1]), grid on
        
        
    end
    
    taudlist(ratccount) = tauD;
end




