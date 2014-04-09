function [M2A M6A M2B M6B] = MSE_IC(aalpha,ggamma,Nsim,K)
% Calculate E[(\hat{\theta}_m(k)-\theta_m)^2] by simulation
% m>0

ttheta_m2 = (aalpha-ggamma)*(aalpha^(2-1));
ttheta_m6 = (aalpha-ggamma)*(aalpha^(6-1));

AIC  = zeros(K+1,1);
BIC  = zeros(K+1,1);
FIC2 = zeros(K+1,1);
FIC6 = zeros(K+1,1);
irf  = zeros(Nsim,4);

T = 200;

model = arima('Constant',0,'AR',aalpha,'MA',ggamma,'Variance',1);

[Y,~,~] = simulate(model,T,'numPaths',Nsim);

counter = 0;
for i=1:Nsim
    
    %Criterion Selection
    
    for k=1:13
        [AIC(k) BIC(k)] = AICBIC(k-1,Y(:,i));
        [FIC2(k)] = FIC(2,k-1,K,Y(:,i));
        [FIC6(k)] = FIC(6,k-1,K,Y(:,i));
    end
    
    [~,orderAIC]  = max(AIC);
    [~,orderBIC]  = max(BIC);
    [~,orderFIC2] = max(FIC2);
    [~,orderFIC6] = max(FIC6);
    orderAIC = orderAIC-1;
    orderBIC = orderBIC-1;
    orderFIC2= orderFIC2-1;
    orderFIC6= orderFIC6-1;
    ICOrder       = [orderAIC orderBIC orderFIC2 orderAIC orderBIC ...
        orderFIC6];
    
    m = [2 2 2 6 6 6];
    
    for s=1:6
        I_k     = [eye(ICOrder(s)-1) zeros(ICOrder(s)-1,1)];
        X       = lagY(Y(:,i),ICOrder(s));
        OLS     = (X'*X)\(X'*Y(ICOrder(s)+1:T,i));
        pphi    = [OLS';I_k];
        irfhold = pphi^m(s);
        irf(i,s)= irfhold(1,1);
    end
    counter = counter+1;
    if(mod(counter,100)==0)
        fprintf('Simulation = %d',counter);
        fprintf('\n')
    end
end

M2A = sqrt(mean((irf(:,3)-ttheta_m2).^2)/mean((irf(:,1)-ttheta_m2).^2));
M2B = sqrt(mean((irf(:,3)-ttheta_m2).^2)/mean((irf(:,2)-ttheta_m2).^2));
M6A = sqrt(mean((irf(:,6)-ttheta_m6).^2)/mean((irf(:,4)-ttheta_m6).^2));
M6B = sqrt(mean((irf(:,6)-ttheta_m6).^2)/mean((irf(:,5)-ttheta_m6).^2));
