% forcingco2.m - CO2-only radiative forcing: F = 5.35*ln(C/C_0)
% Used in: Seshadri (2017) doi:10.1007/s00382-016-3519-3
%

function F = forcingco2(tlist)

global x2_PI tco2 x2 

F = NaN(size(tlist));
for i = 1:numel(tlist)
    t = tlist(i);
    
    x2_t = interp1(tco2,x2,t);
    
    % radiative forcing at time t
    F2_t = 5.35*log(x2_t/x2_PI);
    
    F(i) = F2_t;
    
end


