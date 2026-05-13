% Fig. 4: Contour plot of tolerance theta vs alpha=T/tau and x=t/T
% Exact vs 1st-3rd order series approx. Boundaries are hyperbola families.
%
% Paper: Seshadri (2021), Clim. Dyn., doi:10.1007/s00382-021-05739-3
%

clear, close all

pal = {'k','r','b','g','m','c','y'};

beta = 1;


%% gamma = 0.5

gamma = 0.5;

x = linspace(0.01,0.99,1000);

dx = diff(x); dx = dx(1);

% emissions

m = (x <= 0.5).*(beta*x.^gamma) + (x > 0.5).*(beta*(1-x).^gamma);

% repeated integrals

m1 = cumsum(m)*dx;
m2 = cumsum(m1)*dx;

m3 = cumsum(m2)*dx;
m4 = cumsum(m3)*dx;

m5 = cumsum(m4)*dx;
m6 = cumsum(m5)*dx;

m7 = cumsum(m6)*dx;
m8 = cumsum(m7)*dx;

m9 = cumsum(m8)*dx;
m10 = cumsum(m9)*dx;

%% ratio alpha
alphalist = linspace(0.01,2,2000);


%% Approximate ratio of timescales
ratexactmat = NaN(numel(alphalist),numel(x));
ratmato1 = NaN(numel(alphalist),numel(x));
ratmato2 = NaN(numel(alphalist),numel(x));
ratmato3 = NaN(numel(alphalist),numel(x));
ratmato4 = NaN(numel(alphalist),numel(x));
ratmato5 = NaN(numel(alphalist),numel(x));
ratmato6 = NaN(numel(alphalist),numel(x));
ratmato7 = NaN(numel(alphalist),numel(x));
ratmato8 = NaN(numel(alphalist),numel(x));


for i = 1:numel(alphalist)
    
    alpha = alphalist(i);
    
    %% Exact
    CEi = [0 exp(-alpha*x).*cumsum(exp(alpha*x).*m)*dx];
    dCEi = diff(CEi)/dx;
    ratexact = abs(dCEi./CEi(2:numel(CEi)).*m1./m-1);
       
    
    %% Approximate
    
    rato1i = abs((1 - alpha*m1./m )...
        ./(1 - alpha*m2./m1 )-1);

    rato2i = abs((1 - alpha*m1./m + alpha^2*m2./m )...
        ./(1 - alpha*m2./m1 + alpha^2*m3./m1 )-1);
    
    rato3i = abs((1 - alpha*m1./m + alpha^2*m2./m - alpha^3*m3./m  )...
        ./(1 - alpha*m2./m1 + alpha^2*m3./m1 - alpha^3*m4./m1 )-1);
    
    rato4i = abs((1 - alpha*m1./m + alpha^2*m2./m - alpha^3*m3./m + alpha^4*m4./m )...
        ./(1 - alpha*m2./m1 + alpha^2*m3./m1 - alpha^3*m4./m1 + alpha^4*m5./m1)-1);
    
    rato5i = abs((1 - alpha*m1./m + alpha^2*m2./m - alpha^3*m3./m + alpha^4*m4./m - alpha^5*m5./m )...
        ./(1 - alpha*m2./m1 + alpha^2*m3./m1 - alpha^3*m4./m1 + alpha^4*m5./m1 - alpha^5*m6./m1 )-1);
    
    rato6i = abs((1 - alpha*m1./m + alpha^2*m2./m - alpha^3*m3./m + alpha^4*m4./m - alpha^5*m5./m + alpha^6*m6./m )...
        ./(1 - alpha*m2./m1 + alpha^2*m3./m1 - alpha^3*m4./m1 + alpha^4*m5./m1 - alpha^5*m6./m1 + alpha^6*m7./m1 )-1);

    ratexactmat(i,:) = ratexact';
    
    ratmato1(i,:) = rato1i';
    
    ratmato2(i,:) = rato2i';
    ratmato3(i,:) = rato3i';
    ratmato4(i,:) = rato4i';
    ratmato5(i,:) = rato5i';
    ratmato6(i,:) = rato6i';
 
end


%% Plot figure

% plot different linestyles, for legend
v0 = [0.01];

figure(2), subplot(1,2,1),

contour(alphalist,x,ratexactmat',v0,'k','LineWidth',2,'ShowText','off'); xlabel('dimensionless parameter $\alpha=T / \tau$','Interpreter','latex'), ylabel('rescaled time, $x$','Interpreter','latex'), hold on

contour(alphalist,x,ratmato1',v0,'r+','ShowText','off'), 

contour(alphalist,x,ratmato2',v0,'g.','ShowText','off'), 

contour(alphalist,x,ratmato3',v0,'b','ShowText','off'),

% plot all
v = [0.01 0.05 0.1 0.25];

[c,h] = contour(alphalist,x,ratexactmat',v,'k','LineWidth',2,'ShowText','off'); 
clabel(c,h)

contour(alphalist,x,ratmato1',v,'r+','ShowText','off'), 

contour(alphalist,x,ratmato2',v,'g.','ShowText','off'), 

contour(alphalist,x,ratmato3',v,'b','ShowText','off'),

xlim([0 2]), grid on
legend('Exact','1^{st}-order approx.','2^{nd}-order approx.','3^{rd}-order approx.','Location','NorthEast')
legend('boxoff')

title('$\gamma=0.5$','Interpreter','latex')

%% gamma = 2

gamma = 2;

x = linspace(0.01,0.99,1000);

dx = diff(x); dx = dx(1);

% emissions

m = (x <= 0.5).*(beta*x.^gamma) + (x > 0.5).*(beta*(1-x).^gamma);

% repeated integrals

m1 = cumsum(m)*dx;
m2 = cumsum(m1)*dx;

m3 = cumsum(m2)*dx;
m4 = cumsum(m3)*dx;

m5 = cumsum(m4)*dx;
m6 = cumsum(m5)*dx;

m7 = cumsum(m6)*dx;
m8 = cumsum(m7)*dx;

m9 = cumsum(m8)*dx;
m10 = cumsum(m9)*dx;

%% ratio alpha
alphalist = linspace(0.01,2,2000);


%% Approximate ratio of timescales
ratexactmat = NaN(numel(alphalist),numel(x));
ratmato1 = NaN(numel(alphalist),numel(x));
ratmato2 = NaN(numel(alphalist),numel(x));
ratmato3 = NaN(numel(alphalist),numel(x));
ratmato4 = NaN(numel(alphalist),numel(x));
ratmato5 = NaN(numel(alphalist),numel(x));
ratmato6 = NaN(numel(alphalist),numel(x));
ratmato7 = NaN(numel(alphalist),numel(x));
ratmato8 = NaN(numel(alphalist),numel(x));


for i = 1:numel(alphalist)
    
    alpha = alphalist(i);
    
    %% Exact
    CEi = [0 exp(-alpha*x).*cumsum(exp(alpha*x).*m)*dx];
    dCEi = diff(CEi)/dx;
    ratexact = abs(dCEi./CEi(2:numel(CEi)).*m1./m-1);
       
    
    %% Approximate
    
    rato1i = abs((1 - alpha*m1./m )...
        ./(1 - alpha*m2./m1 )-1);

    rato2i = abs((1 - alpha*m1./m + alpha^2*m2./m )...
        ./(1 - alpha*m2./m1 + alpha^2*m3./m1 )-1);
    
    rato3i = abs((1 - alpha*m1./m + alpha^2*m2./m - alpha^3*m3./m  )...
        ./(1 - alpha*m2./m1 + alpha^2*m3./m1 - alpha^3*m4./m1 )-1);
    
    rato4i = abs((1 - alpha*m1./m + alpha^2*m2./m - alpha^3*m3./m + alpha^4*m4./m )...
        ./(1 - alpha*m2./m1 + alpha^2*m3./m1 - alpha^3*m4./m1 + alpha^4*m5./m1)-1);
    
    rato5i = abs((1 - alpha*m1./m + alpha^2*m2./m - alpha^3*m3./m + alpha^4*m4./m - alpha^5*m5./m )...
        ./(1 - alpha*m2./m1 + alpha^2*m3./m1 - alpha^3*m4./m1 + alpha^4*m5./m1 - alpha^5*m6./m1 )-1);
    
    rato6i = abs((1 - alpha*m1./m + alpha^2*m2./m - alpha^3*m3./m + alpha^4*m4./m - alpha^5*m5./m + alpha^6*m6./m )...
        ./(1 - alpha*m2./m1 + alpha^2*m3./m1 - alpha^3*m4./m1 + alpha^4*m5./m1 - alpha^5*m6./m1 + alpha^6*m7./m1 )-1);

    ratexactmat(i,:) = ratexact';
    
    ratmato1(i,:) = rato1i';
    
    ratmato2(i,:) = rato2i';
    ratmato3(i,:) = rato3i';
    ratmato4(i,:) = rato4i';
    ratmato5(i,:) = rato5i';
    ratmato6(i,:) = rato6i';
 
end


%% Plot figure

% plot different linestyles, for legend
v0 = [0.01];

figure(2), subplot(1,2,2),

contour(alphalist,x,ratexactmat',v0,'k','LineWidth',2,'ShowText','off'); xlabel('dimensionless parameter $\alpha=T / \tau$','Interpreter','latex'), ylabel('rescaled time, $x$','Interpreter','latex'), hold on

contour(alphalist,x,ratmato1',v0,'r+','ShowText','off'), 

contour(alphalist,x,ratmato2',v0,'g.','ShowText','off'), 

contour(alphalist,x,ratmato3',v0,'b','ShowText','off'),

% plot all
v = [0.01 0.05 0.1 0.25];

[c,h] = contour(alphalist,x,ratexactmat',v,'k','LineWidth',2,'ShowText','off'); 
clabel(c,h)

contour(alphalist,x,ratmato1',v,'r+','ShowText','off'), 

contour(alphalist,x,ratmato2',v,'g.','ShowText','off'), 

contour(alphalist,x,ratmato3',v,'b','ShowText','off'),

xlim([0 2]), grid on
legend('Exact','1^{st}-order approx.','2^{nd}-order approx.','3^{rd}-order approx.','Location','NorthEast')
legend('boxoff')

title('$\gamma=2.0$','Interpreter','latex')

