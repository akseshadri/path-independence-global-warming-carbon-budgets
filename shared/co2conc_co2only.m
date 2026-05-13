% co2conc_co2only.m - CO2 concentration from emissions via Joos et al. (2013)
%   h(t) = 0.276*exp(-t/4.3) + 0.282*exp(-t/36.5) + 0.224*exp(-t/394) +
%   0.218
% Inputs: gr (growth rate), tstab (stabilization time), taum2 (decay timescale)
% Outputs: timp (time), xco2 (concentration ppm), m1 (cumulative emissions)
% Ref: Joos et al. (2013), doi:10.5194/acp-13-2793-2013
% Used in: Seshadri (2017) and Seshadri (2021)
%

function [timp,xco2,m1] = co2conc_co2only(gr,tstab,taum2)

global mu1 mu2 mu3 mu4 tau1 tau2 tau3

[tpast,mpast] = textread('co2.txt');

mpast = mpast*3.67*1e9; % kg CO2
mpast = mpast*120/9.5e14; % ppmv CO2

tpast = tpast - max(tpast);

tfuture = linspace(1,1000,1000)';
m0 = mpast(numel(mpast));

gdp = (1 + gr).^min(tfuture,tstab);
mfuture = m0*gdp.*exp(-tfuture/taum2);

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
