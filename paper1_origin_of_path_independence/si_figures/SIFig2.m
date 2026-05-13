% SI Fig. 2: Spread of damping timescale tau_D across 16 CMIP5 GCMs
% Based on Geoffroy et al. (2013b) EBM parameter estimates
%
% Paper: Seshadri (2017), Clim. Dyn. 49:3383-3401, doi:10.1007/s00382-016-3519-3
%


clear, close all

global RFflag cs cd beta nu gamma 

RFflag = 2;

palexact = {'k','b','r'};
palapprox = {'k--','b--','r--'};

modellist= {'BCC-CSM1.1','BNU-ESM','CanESM2','CCSM4','CNRM-CM5.1','CSIRO-Mk3.6.0','FGOALS-s2','GFDL-ESM2M','GISS-E2-R','HadGEM2-ES','INM-CM4','IPSL-CM5A-LR'...
    ,'MIROC5','MPI-ESM-LR','MRI-CGCM3','NorESM1-M'}; 

%% read heat capacities and gamma
dat1 = textread('heatcapacity.txt');
cslist = dat1(:,1); cdlist = dat1(:,2); gammalist = dat1(:,3);

%% read other parameters
dat2 = textread('betaandnu.txt');
betalist = dat2(:,2); nulist = dat2(:,3);



%% Simulation loop
memat = NaN(16,1);
mtmat = NaN(16,1);

tauDlist = NaN(16,1);
tau1list = NaN(16,1);

for n = 1:16
    
    
    %% Assign parameters
    cs = cslist(n);
    cd = cdlist(n);
    epsilon = cs/cd;
    r = 1/epsilon;
    beta = betalist(n);
    gamma = gammalist(n);
    nu = nulist(n);
      
    %% taud and tau1, tau2
    tauD = cd*(beta+nu*gamma)/nu/gamma^2;
    tau1 = cd*(beta+nu*gamma)/beta/gamma;
    tau2 = cs/(beta+nu*gamma);
    
    tauDlist(n) = tauD;
    tau1list(n) = tau1; 
    
    
end

%% plot fig 4
figure(2), 

i = 1:16 ~= 11;


[n,td] = hist(tauDlist(i),4);
subplot(1,2,2),
plot(td,n, 'k+-'), xlabel('damping-timescale \tau_D   (years)'), ylabel('histogram count'),
    subplot(1,2,1), plot(tauDlist, 'k+'), xlabel('GCM'), ylabel('damping-timescale \tau_D (years)'),...
        xlim([1 16]), ylim([100 1600])
  
