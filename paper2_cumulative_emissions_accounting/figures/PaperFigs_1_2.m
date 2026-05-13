% Figs. 1-2: Path independence conditions for CO2 in EBM
% Fig 1: (a) Emissions (b) Warming vs cum. emissions (c) Directional derivatives (d) Timescales
% Fig 2: Relative effects bounded by ratio of timescales (Cf. Eq. 13)
%
% Paper: Seshadri (2021), Clim. Dyn., doi:10.1007/s00382-021-05739-3
%

clear, close all

global RFflag gr cs cd beta nu gamma tco2 x2 tstab mu1 mu2 mu3 mu4 tau1 tau2 tau3 x2_PI

RFflag = 2;

pal = {'k','b','r','c','g','m'};
palpoint = {'.','.','.','.','.','.'};

taum2list = [20 inf 100 10 50 100];
ratc = 1/0.075;

%% conversion factor
cfp = 1/(3.67*1e9*120/9.5e14)/1e3;

%% economy parameters
grlist = [0.02 0.025 0.03 0.02 0.035 0.025];
tstablist = [50 50 50 100 100 100];

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
cd0 = cs0*ratc;

cs = cs0; cd = cd0;

gamma = 0.67; % W/m^2.K
nu = 1.28; % efficacy of ocean heat uptake
DeltaTECS = 3.1; % K

RFnu = 5.35; % W/m^2
beta = RFnu*log(2)/DeltaTECS; % W/m^2.K

Ts0 = 0.0; Td0 = 0.0;

%% taud and taus, tauf
tauD = cd*(beta+nu*gamma)/nu/gamma^2; % years
taus = cd*(beta+nu*gamma)/beta/gamma; % years
tauf = cs/(beta+nu*gamma); % years

%% Conversion from u to Ts
conv_u_to_Ts = RFnu*tauf/cs; % conversion factor from dimensionless u to Ts in K

% required for computing gradients

%% ODE options
optionsc = odeset('RelTol',1e-2,'AbsTol',1e5,'MaxStep',1.0);
optionst = odeset('RelTol',1e-3,'AbsTol',[1e-4 1e-4],'MaxStep',0.2);


for maincount = 1:6
    
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
    
    % interpolate to tco2
    Tsp = interp1(s,Ts,tco2);
    
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
    
    %% Calculate m1/g, r/g and g/g0
    m1byg = m1t1./gt1;
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
        
        % timescales at t
        taum1i = taum1(i);
        tauri = taur(i);
        
        % variables at t
        ri = rt1(i);
        gi = gt1(i);
        m1i = m1t1(i);
        
        % ratios 
        m1bygi = m1byg(1:i);
        rbygi = rbyg(1:i);
        gbyg0i = gbyg0(1:i);
        
        % time variables
        zi = t1(1:i);
        ti = t1(i);
        dz = t1(2) - t1(1);
        
        %% gradient
        % derivative w.r.t. t
        delu_by_delt = (1/tauD)*(log(gbyg0(i)) - 1/taus*exp(-ti/taus)*sum(exp(zi/taus).*log(gbyg0i)*dz));

        % derivative w.r.t r
        delu_by_delr = m1i/gi+1/tauD*exp(-ti/taus)*sum(exp(zi/taus).*m1bygi*dz);
        
        % derivative w.r.t. m1
        delu_by_delm1 =  ri/gi+1/tauD*exp(-ti/taus)*sum(exp(zi/taus).*rbygi*dz);
        
        %% direction vector
        % change in r
        dr_by_dt = abs(ri./tauri);
        
        % change in m1
        dm1_by_dt = abs(m1i/taum1i);

        %% directional derivative
        
        % t
        dd_t(i) = conv_u_to_Ts*delu_by_delt; 
        % r
        dd_r(i) = conv_u_to_Ts*delu_by_delr*dr_by_dt;
        % m1
        dd_m1(i) = conv_u_to_Ts*delu_by_delm1*dm1_by_dt;
        
    end
    
    
    %% Plot figure 1
    
    figure(1),
    
    % emissions vs time
    subplot(2,2,1),
    plot(t1,m*cfp, char(pal(maincount)), 'LineWidth', 1.0), hold on, ...
        xlabel('year','Interpreter','latex'), ylabel('emissions $(PgC.year^{-1})$','Interpreter','latex'),...
        xlim([1900 2200]), grid on
    
    % global warming vs. cumulative co2 emissions
    subplot(2,2,2),
    plot(m1*cfp,Tsp, char(pal(maincount)), 'LineWidth', 1.0), hold on, ...
        ylabel('global warming $(K)$','Interpreter','latex'), xlabel('cumulative emissions $(PgC)$','Interpreter','latex'),...
       xlim([0 4000]), ylim([0 5]), grid on
    
    % directional derivative
    pal1 = strcat(char(pal(maincount)),'--');
    pal2 = strcat(char(pal(maincount)),':');
    pal3 = strcat(char(pal(maincount)),'-');
    
    subplot(2,2,3), semilogy(t1,dd_t,pal1, 'LineWidth',0.75), hold on, ...
        semilogy(t1,dd_r,pal2, 'LineWidth',0.75),...
        semilogy(t1,dd_m1,pal3, 'LineWidth',1.0)
        xlabel('year','Interpreter','latex'), ylabel('component of directional derivative $(year^{-1})$','Interpreter','latex'), ...
        xlim([1900 2200]), ylim([1e-5 1e-1]), grid on
    yticks([1e-5 1e-4 1e-3 1e-2 1e-1])
    
    legend('time','airborne fraction','cumulative emissions','Location','SouthEast')
    legend('boxoff')
    
    % time-scale for evolution
    
    taur = abs(taur);
    subplot(2,2,4), 
    semilogy(t1,taur,pal2,'LineWidth',0.75), hold on,...
        semilogy(t1,taum1,pal3,'LineWidth',1.0),...
    xlabel('year','Interpreter','latex'), ylabel('timescale for evolution (years)','Interpreter','latex'), ...
    xlim([1900 2200]), ylim([10 10^4])  
    legend('airborne fraction','cumulative emissions','Location','SouthEast')
    legend('boxoff'), grid on
    
    
    %% Plot Figure 2
    
    % airborne fraction / cumulative emissions
    rat_dd_r_by_m1 = (dd_r)./dd_m1;
    rat_taum1_by_taur = taum1./taur;
    
    % time / cumulative emissions
    rat_dd_t_by_m1 = (dd_t)./dd_m1;
    rat_taum1_by_tauD = taum1./tauD;
    
    pal4 = strcat(char(pal(maincount)),char(palpoint(maincount)));
    
    xlin = [1e-3 1];
    ylin = xlin;
    
    figure(2), 
    
    subplot(1,2,1), loglog(rat_taum1_by_taur,rat_dd_r_by_m1, pal4,'LineWidth',0.5,'MarkerSize',3), hold on,...
        xlim([1e-3 1]), ylim([1e-3 1]), xlabel('ratio of timescales $\left|\tau_{M}/\tau_{r}\right|$','Interpreter','latex'),...
        ylabel('relative effect','Interpreter','latex'), title('airborne fraction / cumulative emissions','Interpreter','latex'),...
        plot(xlin,ylin,'k--','LineWidth',1), grid on,...
        xlim([1e-3 1]), ylim([1e-3 1]),
    
        
    subplot(1,2,2), loglog(rat_taum1_by_tauD,rat_dd_t_by_m1, pal4,'LineWidth',0.5,'MarkerSize',3), hold on,...
        xlim([1e-3 1]), ylim([1e-3 1]), xlabel('ratio of timescales $\tau_M / \tau_D$','Interpreter','latex'),...
        title('time'), title('time / cumulative emissions','Interpreter','latex'),...
        plot(xlin,ylin,'k--','LineWidth',1), grid on,...
        xlim([1e-3 1]), ylim([1e-3 1]),
    
end







