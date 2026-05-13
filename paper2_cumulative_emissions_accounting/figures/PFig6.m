% Fig. 6: HFC143a global warming vs cumulative emissions (tau=52 yr)
% Left: slow mitigation (path independence fails) Right: rapid mitigation (holds)
%
% Paper: Seshadri (2021), Clim. Dyn., doi:10.1007/s00382-021-05739-3
%

clear, close all

global cs cd beta nu gamma years rf

pal = {'g','b','r','c','k','m'};
pal2 = {'k','k--','k.-','k:','k+-'};

%% 2-box model parameters
cs0 = 8.2; % W*yr / m^2.K
ratc = 1/0.075;
cd0 = cs0*ratc;

cs = cs0; cd = cd0;

gamma = 0.67; % W/m^2.K
nu = 1.28; % efficacy of ocean heat uptake
DeltaTECS = 3.1; % K

RFnu = 5.35; % W/m^2
beta = RFnu*log(2)/DeltaTECS; % W/m^2.K

Ts0 = 0.0; Td0 = 0.0;

%% read past HFC emissions
hfcpast = xlsread('hfchist.xlsx');

yearsp = hfcpast(:,1);
hfcp1 = hfcpast(:,2);
hfcp2 = hfcpast(:,3);

%% read RCP scenarios
rcp3 = xlsread('rcp3proj.xlsx');
rcp45 = xlsread('rcp45proj.xlsx');
rcp6 = xlsread('rcp6proj.xlsx');
rcp85 = xlsread('rcp85proj.xlsx');

yearsf3 = rcp3(:,1);
hfcf3_1 = rcp3(:,2);
hfcf3_2 = rcp3(:,3);

yearsf45 = rcp45(:,1);
hfcf45_1 = rcp45(:,2);
hfcf45_2 = rcp45(:,3);

yearsf6 = rcp6(:,1);
hfcf6_1 = rcp6(:,2);
hfcf6_2 = rcp6(:,3);

yearsf85 = rcp85(:,1);
hfcf85_1 = rcp85(:,2);
hfcf85_2 = rcp85(:,3);

%% Interpolate
years = 1991:2200;

% RCP 3

ygiven_3 = [yearsp(1:25)' yearsf3(10:20)'];
hfc1given_3 = [hfcp1(1:25)' hfcf3_1(10:20)'];
hfc2given_3 = [hfcp2(1:25)' hfcf3_2(10:20)'];

hfc1_3 = interp1(ygiven_3,hfc1given_3,years);
hfc2_3 = interp1(ygiven_3,hfc2given_3,years);

% RCP 4.5

ygiven_45 = [yearsp(1:25)' yearsf45(10:49)'];
hfc1given_45 = [hfcp1(1:25)' hfcf45_1(10:49)'];
hfc2given_45 = [hfcp2(1:25)' hfcf45_2(10:49)'];

hfc1_45 = interp1(ygiven_45,hfc1given_45,years);
hfc2_45 = interp1(ygiven_45,hfc2given_45,years);

% RCP 6

ygiven_6 = [yearsp(1:25)' yearsf6(10:51)'];
hfc1given_6 = [hfcp1(1:25)' hfcf6_1(10:51)'];
hfc2given_6 = [hfcp2(1:25)' hfcf6_2(10:51)'];

hfc1_6 = interp1(ygiven_6,hfc1given_6,years);
hfc2_6 = interp1(ygiven_6,hfc2given_6,years);

% RCP 8.5

ygiven_85 = [yearsp(1:25)' yearsf85(10:31)'];
hfc1given_85 = [hfcp1(1:25)' hfcf85_1(10:31)'];
hfc2given_85 = [hfcp2(1:25)' hfcf85_2(10:31)'];

hfc1_85 = interp1(ygiven_85,hfc1given_85,years);
hfc2_85 = interp1(ygiven_85,hfc2given_85,years);

%% Only RCP 3 and 8.5 scenarios are smooth, use them

hfc1_4 = 0.75*hfc1_3 + 0.25*hfc1_85;
hfc2_4 = 0.75*hfc2_3 + 0.25*hfc2_85;

hfc1_6 = 0.25*hfc1_3 + 0.75*hfc1_85;
hfc2_6 = 0.25*hfc2_3 + 0.75*hfc2_85;

%% Molecular weights and lifetimes
Matm = 5.15e18; % kg
Mair = 28.97; % g/mol
Mhfc2 = 84.04; % g/mol
tau = 52;


%% Emissions scenarios
Alist = [0 60 -60 30 -30];
taulist = [400 400 400 400 400];
clist = [30 30 30 30 30];
%clist = [5 10 0 0 0];

emismat = NaN(4,210);

for i = 1:5
    Ai = Alist(i);
    taui = taulist(i);
    ci = clist(i);
    yi = 1:210; yi = max(0,yi-ci);
    emisi =  hfc2_3 + Ai*sin(2*pi*yi/taui);
    emismat(i,:) = emisi;
end

%% Calculate global warming

for i = 1:4
    
    % emis scenario
    emis = emismat(i,:); % Ggrams
   
    % concentration
    [tout,conc1,cumemis1] = getconc(emis,tau); % concentration is in Ggrams = 1e6 kg
    
    concppbv = conc1*1e6/Matm*Mair/Mhfc2*1e9; % in ppbv
    
    % radiative forcing
    rf = concppbv*0.13; % W/m^2
    
    % global warming
    optionsc = odeset('RelTol',1e-2,'AbsTol',1e5,'MaxStep',1.0);
    optionst = odeset('RelTol',1e-3,'AbsTol',[1e-4 1e-4],'MaxStep',0.2);
    
    mintime = 1991; maxtime = 2200;
    Ts0 = 0.0; Td0 = 0.0;
    
    [s,Y] = ode45(@twoboxmodelsingleforcer,[mintime maxtime],[Ts0 Td0],optionst);
    
    Ts = Y(:,1)-Ts0; [Tsmax,indexmax] = max(Ts);
    
    Td = Y(:,2)-Td0;
    
    Tsy = interp1(s,Ts,years);
    M1 = cumemis1/1e6; % Gtonnes or Petagrams

    % plot
    figure(5),
    subplot(2,2,1), plot(years,emis,char(pal(i)),'LineWidth',1.5), hold on,...
        xlabel('year','Interpreter','latex'), ylabel('HFC143a emis. (Gg / year)','Interpreter','latex'),...
        xlim([1990 2100]), ylim([0 150]), grid on
    subplot(2,2,3), plot(M1,Tsy,char(pal(i)),'LineWidth',1.5), hold on,...
        xlabel('cumulative emissions (Pg)','Interpreter','latex'), ylabel('global warming (K)','Interpreter','latex'),...
        grid on
end

%% Alternative scenario family
taulistmit = [15 12 9 6 3];
clist = [30 37 44 51 58];

for i = 1:4
    
    % emis scenario
    ci = clist(i);
    taumiti = taulistmit(i);
    yi = 1:210; yi = max(0,yi-ci);
    fi = exp(-yi/taumiti);
    
    emisi = emismat(i,:);
    emis = emisi.*fi; % Ggrams
   
    % concentration
    [tout,conc1,cumemis1] = getconc(emis,tau); % concentration is in Ggrams = 1e6 kg
    
    concppbv = conc1*1e6/Matm*Mair/Mhfc2*1e9; % in ppbv
    
    % radiative forcing
    rf = concppbv*0.13; % W/m^2
    
    % global warming
    optionsc = odeset('RelTol',1e-2,'AbsTol',1e5,'MaxStep',1.0);
    optionst = odeset('RelTol',1e-3,'AbsTol',[1e-4 1e-4],'MaxStep',0.2);
    
    mintime = 1991; maxtime = 2200;
    Ts0 = 0.0; Td0 = 0.0;
    
    [s,Y] = ode45(@twoboxmodelsingleforcer,[mintime maxtime],[Ts0 Td0],optionst);
    
    Ts = Y(:,1)-Ts0; [Tsmax,indexmax] = max(Ts);
    
    Td = Y(:,2)-Td0;
    
    Tsy = interp1(s,Ts,years);
    M1 = cumemis1/1e6; % Gtonnes or Petagrams

    % plot
    figure(5),
    subplot(2,2,2), plot(years,emis,char(pal2(i)),'LineWidth',1.5), hold on,...  
        xlabel('year','Interpreter','latex'),...
        xlim([1990 2100]), ylim([0 150]), grid on
    subplot(2,2,4), plot(M1,Tsy,char(pal2(i)),'LineWidth',1.5), hold on,...
        xlabel('cumulative emissions (Pg)','Interpreter','latex'), ...
        grid on
end


subplot(2,2,1), title('Slow Mitigation','Interpreter','latex'),...
    subplot(2,2,2), title('Rapid Mitigation','Interpreter','latex'),...