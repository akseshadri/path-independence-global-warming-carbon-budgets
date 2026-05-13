% twoboxmodelco2.m - Two-box EBM driven by CO2 forcing only
%
% F(t) = 5.35 * ln(C(t)/C_0) with C_0 = preindustrial CO2.
% Used in: Seshadri (2017) doi:10.1007/s00382-016-3519-3
%          Seshadri (2021) doi:10.1007/s00382-021-05739-3
%

function dy = twoboxmodelco2(t,y)

global x2_PI cs cd beta nu gamma tco2 x2 

% concentration at time t
x2_t = interp1(tco2,x2,t);

% radiative forcing at time t
F = 5.35*log(x2_t/x2_PI);

dy = zeros(2,1);    % a column vector
dy(1) = -(beta+nu*gamma)/cs*y(1) + nu*gamma/cs*y(2) + F/cs;
dy(2) = gamma/cd*y(1) - gamma/cd*y(2);
