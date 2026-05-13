% Concentration for single-lifetime species: h(t)=exp(-t/tau). For HFCs, N2O etc.
% Paper: Seshadri (2021), doi:10.1007/s00382-021-05739-3
%

function [tout,conc,cumemis] = getconc(emis,tau)

cumemis = cumsum(emis);

timp = 0:250;

imp = exp(-timp/tau);

timp = timp + 1991; 

%% concentration
conc = conv(emis,imp); 

%% output
tout = timp(1:210);
conc = conc(1:210);
cumemis = cumemis(1:210)';
