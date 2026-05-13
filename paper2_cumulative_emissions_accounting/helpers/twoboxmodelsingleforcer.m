% Two-box EBM for a single forcer with precomputed RF time series
% Paper: Seshadri (2021), doi:10.1007/s00382-021-05739-3
%

function dy = twoboxmodelsingleforcer(t,y)

global cs cd beta nu gamma years rf

F = interp1(years,rf,t);

dy = zeros(2,1);    % a column vector
dy(1) = -(beta+nu*gamma)/cs*y(1) + nu*gamma/cs*y(2) + F/cs;
dy(2) = gamma/cd*y(1) - gamma/cd*y(2);
