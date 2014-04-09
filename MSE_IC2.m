function [M2] = MSE_IC2(aalpha,ggamma,Nsim,K)
% Calculate E[(\hat{\theta}_m(k)-\theta_m)^2] by simulation
% m>0

ttheta_m2 = (aalpha-ggamma)*(aalpha^(2-1));

AIC  = zeros(K+1,1);
BIC  = zeros(K+1,1);
FIC2 = zeros(K+1,1);
irf  = zeros(Nsim,3);

T = 200;

model = arima('Constant',0,'AR',aalpha,'MA',ggamma,'Variance',1);

[Y,~,~] = simulate(model,T,'numPaths',Nsim);

counter = 0;
for i=1:Nsim
    
    %Criterion Selection
    
    for k=1:13
        [AIC(k) BIC(k)] = AICBIC(k-1,Y(:,i));
        [FIC2(k)] = criterion(2,k-1,K,Y(:,i));
    end
    
    [~,orderAIC]  = max(AIC);
    [~,orderBIC]  = max(BIC);
    [~,orderFIC2] = max(FIC2);
    orderAIC = orderAIC-1;
    orderBIC = orderBIC-1;
    orderFIC2= orderFIC2-1;
    ICOrder       = [orderAIC orderBIC orderFIC2];
    
    m=2;
    for s=1:3
        I_k     = [eye(ICOrder(s)-1) zeros(ICOrder(s)-1,1)];
        X       = lagY(Y(:,i),ICOrder(s));
        OLS     = (X'*X)\(X'*Y(ICOrder(s)+1:T,i));
        pphi    = [OLS';I_k];
        irfhold = pphi^m;
        irf(i,s)= irfhold(1,1);
    end
    counter = counter+1;
    if(mod(counter,100)==0)
        fprintf('Simulation = %d',counter);
        fprintf('\n')
    end
end

M2 = sqrt(mean((irf(:,3)-ttheta_m2).^2)/mean((irf(:,1)-ttheta_m2).^2));