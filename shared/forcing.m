% forcing.m - Compute total radiative forcing at specified times
% Returns F(t) = F_CO2 + F_other. Historical for t<=2014, projections after.
% Used in: Seshadri (2017) doi:10.1007/s00382-016-3519-3
%

function F = forcing(tlist)

global RFflag F01 x01 x2_PI tbc x1 tco2 x2 thist Fhist tpd Fpd t45 F45

F = NaN(size(tlist));
for i = 1:numel(tlist)
    t = tlist(i);
    % concentration at time t
    if t <= 2014
        F(i) = interp1(thist,Fhist,t);
    else
        % concentration at time t
        x1_t = interp1(tbc,x1,t);
        x2_t = interp1(tco2,x2,t);
        
        % radiative forcing at time t
        F2_t = 5.35*log(x2_t/x2_PI);
        F1_t = F01*x1_t/x01;
        
        if RFflag == 1
            Fothers_t = interp1(tpd,Fpd,t);
        else
            Fothers_t = interp1(t45,F45,t);
        end
        
        F(i) = F1_t + F2_t + Fothers_t;
    end
end
