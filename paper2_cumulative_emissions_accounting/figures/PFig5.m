% Fig. 5: |tau_M/tau_r| at emissions peak vs alpha=T/tau
% (a) Normalized shape functions (b) Ratio: markers=exact, lines=cubic approx.
%
% Paper: Seshadri (2021), Clim. Dyn., doi:10.1007/s00382-021-05739-3
%

clear, close all

% palette etc.
pal = {'r','k','b'};
palpoint = {'ro','kx','b+'};

%% Calculate solution from cubic polynomial

% theta
thetalist = logspace(-3,0,20);

% gamma
gammalist = [0.5 1.0 2.0];


for i = 1:numel(gammalist)
    gi = gammalist(i);
    alphastar = NaN(numel(thetalist),1);
    for j = 1:numel(thetalist)
        thetai = thetalist(j);
        
        a3i = 3+thetai*(gi+1);
        a2i = -(2+thetai*(gi+1))*(gi+4);
        a1i = (1+thetai*(gi+1))*(gi+3)*(gi+4);
        a0i = -(gi+1)*(gi+2)*(gi+3)*(gi+4)*thetai;
        
        
        pi = [a3i a2i a1i a0i];
        ri = roots(pi);
        for k = 1:3
            rik = ri(k);
            if isreal(rik)
                ystar = rik; % real root of cubic
            end
        end
        
        alphastar(j) = 2*ystar; % alphastar = ystar/x; x = 1/2
        
        
    end
    % plot
    figure(3),
    subplot(1,2,2), loglog(alphastar,thetalist,char(pal(i)),'LineWidth',1.5), hold on, 
end

%% Calculate ratio from integrals
clear all

% palette etc.
pal = {'r','k','b'};
palpoint = {'ro','kx','b+'};

% parameter ranges
Tlist = 20:100:220;
gammalist = [0.5 1 2.0];
taulist = logspace(1,3,6);

f1peak = NaN(numel(Tlist),numel(gammalist),numel(taulist));

% Loop
x = linspace(0,1,1001);
for i = 1:numel(Tlist)
    T = Tlist(i);
    for j = 1:numel(gammalist)
        gamma = gammalist(j);
        
        t = T*x;
        dt = t(2) - t(1);
        
        % emissions
        beta = 1/0.5^gamma;
        
        m = beta*x.^gamma.*(x<0.5)+beta*(1-x).^gamma.*(x>=0.5);
        
        % cum. emissions
        M = cumsum(m)*dt;
        
        % cum. emissions timescale
        tauM = M./m;
        
        for k = 1:numel(taulist)
            tau = taulist(k);
            
            % excess conc.
            Ce = exp(-t/tau).*cumsum(exp(t/tau).*m)*dt;
            
            % airborne fraction of cumulative emissions
            r = Ce./M;
            
            % airborne fraction timescale
            dr = diff(r)./diff(t);
            r2 = getnpointavg(r,2);
            t2 = getnpointavg(t,2);
            
            taur = abs(r2./dr);
            
            % ratio of timescales at the peak emissions
            tauM2 = getnpointavg(tauM,2);
            f1 = tauM2./taur;
            
            f1peak(i,j,k) = (f1(500)+f1(501))/2;
            
            % Plot
            
            figure(3),
            subplot(1,2,2),
            xplot = T/tau;
            yplot = f1peak(i,j,k);
            
            loglog(xplot,yplot,char(palpoint(j)),'MarkerSize',6), hold on, grid on
            xlim([1e-2 1e1]), ylim([1e-3 1])
            
        end
        
       
        
        % Plot
        
        figure(3),
        subplot(1,2,1), plot(t,m,char(pal(j)),'LineWidth',1.5), hold on,
            
    
    
    end
end



%% Horizontal line and Figure info.

figure(3), subplot(1,2,1), xlabel('years from start of emission','Interpreter','latex'), ylabel('emissions (normalized units)','Interpreter','latex')
subplot(1,2,2), xlabel('dimensionless parameter $\alpha = T/\tau$','Interpreter','latex'), ylabel('ratio of timescales $\left|\tau_{M}/\tau_{r}\right|$','Interpreter','latex')


xline = [1e-2 1e2]; yline = [1e-1 1e-1];

subplot(1,2,2), plot(xline,yline,'k--','LineWidth',1.0), ...
leg = strcat('\gamma=',num2str(gammalist',3)); legend(leg,'Location','NorthWest'), legend('boxoff')
