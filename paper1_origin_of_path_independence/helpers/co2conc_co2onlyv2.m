% CO2 concentration variant: piecewise-linear emissions (peak & decline)
% Uses Joos et al. (2013) impulse response. See shared/co2conc_co2only.m for base version.
% Paper: Seshadri (2017), doi:10.1007/s00382-016-3519-3
%

function [timp,xco2,m,m1] = co2conc_co2onlyv2(Tp,Tf,Mco2,sigmamco2)

% prepare
[tpast,mpast] = textread('co2.txt');

mpast = mpast*3.67*1e9; % kg CO2
mpast = mpast*120/9.5e14; % ppmv CO2

tpast = tpast - max(tpast);

tfuture = linspace(1,1000,1000)';
m0 = mpast(numel(mpast));

% calc
mp = (2*Mco2-m0*Tp)/Tf*(1+sigmamco2*rand(1));

mfuture = (m0 + (mp-m0)/Tp*tfuture).*(tfuture <= Tp) + (mp-(mp-0)*(tfuture-Tp)/(Tf-Tp)).*(tfuture > Tp & tfuture <= Tf) + 0*(tfuture >= Tf);

t = [tpast; tfuture];

m = [mpast; mfuture];
m1 = cumsum(m);

timp = linspace(0,1500,1501);

imp = 0.276*exp(-timp/4.30)+0.282*exp(-timp/36.5)+0.224*exp(-timp/394)+0.217;

timp = timp + 1751; 

%% concentration
xco2 = conv(m,imp) + 285; % ppm

%% output
timp = timp(1:1250);
xco2 = xco2(1:1250);
m1 = m1(1:1250)';
m = m(1:1250);
