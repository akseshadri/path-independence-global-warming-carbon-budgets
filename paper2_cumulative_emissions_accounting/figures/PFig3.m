% Fig. 3: Cubic polynomial g(y) for path independence condition (Eq. 32)
% g(y) is increasing, g(0)<0, single positive root y* in (0,1). Path independence: y < y*.
%
% Paper: Seshadri (2021), Clim. Dyn., doi:10.1007/s00382-021-05739-3
%

clear, close all

% palette
pal1 = {'r','k','b'};
pal2 = {':','-','.-'};

% input parameters

thetalist = [0.05 0.1 0.2];
gammalist = [0.5 1.0 2.0];

y = linspace(-5,5,1000);

count = 0;

for i = 1:2
    theta = thetalist(i);
    for j = 1:3
        
        % calculate
        gamma = gammalist(j);
        c3 = 3+theta*(gamma+1);
        c2 = (2+theta*(gamma+1))*(gamma+4);
        c1 = (1+theta*(gamma+1))*(gamma+4)*(gamma+3);
        c0 = -theta*(gamma+4)*(gamma+3)*(gamma+2)*(gamma+1);
        
        g = c3*y.^3 + c2*y.^2 + c1*y.^1 + c0;
        
        % plot
        pal = strcat(char(pal1(j)),char(pal2(i)));
        figure(3), plot(y,g,pal,'LineWidth',1.5), hold on,
        
    end
end

grid on, 
leg = {'\theta=0.05; \gamma=0.5','\theta=0.05; \gamma=1.0','\theta=0.05; \gamma=2.0',...
    '\theta=0.1; \gamma=0.5','\theta=0.1; \gamma=1.0','\theta=0.1; \gamma=2.0'};%,...
    %'\theta=0.2; \gamma=0.5','\theta=0.2; \gamma=1.0','\theta=0.2; \gamma=2.0'};

legend(leg,'Location','NorthWest')
xlabel('y','Interpreter','latex'),ylabel('g(y)','Interpreter','latex'),
legend('boxoff')
xlim([-4 4]), ylim([-150 150])
        