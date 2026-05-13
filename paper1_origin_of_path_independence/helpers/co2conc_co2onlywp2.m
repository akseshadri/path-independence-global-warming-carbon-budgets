% CO2 concentration variant: emissions with multiplicative 2nd-phase perturbation
% Uses Joos et al. (2013) impulse response. See shared/co2conc_co2only.m for base version.
% Paper: Seshadri (2017), doi:10.1007/s00382-016-3519-3
%

function [timp,xco2,m1] = co2conc_co2onlywp2(gr,tstab,taum2)

global mu1 mu2 mu3 mu4 tau1 tau2 tau3 tn1 tn2 sn

[tpast,mpast] = textread('co2.txt');

mpast = mpast*3.67*1e9; % kg CO2
mpast = mpast*120/9.5e14; % ppmv CO2

tpast = tpast - max(tpast);

tfuture = linspace(1,1000,1000)';
m0 = mpast(numel(mpast));

nc = (tfuture >= tn1 & tfuture <= tn2).*(1+sn*gr).^(tfuture-tn1)+~(tfuture >= tn1 & tfuture <= tn2); 

gdp = (1 + gr).^min(tfuture,tstab);
mfuture = m0*gdp.*exp(-tfuture/taum2).*nc; 

t = [tpast; tfuture];

m = [mpast; mfuture];
m1 = cumsum(m);

timp = linspace(0,1500,1501);

imp = mu1*exp(-timp/tau1)+mu2*exp(-timp/tau2)+mu3*exp(-timp/tau3)+mu4;

timp = timp + 1751; 

%% concentration
xco2 = conv(m,imp) + 285; % ppm

%% output
timp = timp(1:1250);
xco2 = xco2(1:1250);
m1 = m1(1:1250)';