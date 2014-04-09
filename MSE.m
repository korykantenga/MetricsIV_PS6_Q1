function M = MSE(aalpha,ggamma,m,k,Nsim)
% Calculate E[(\hat{\theta}_m(k)-\theta_m)^2] by simulation
% m>0

ttheta_m = (aalpha-ggamma)*(aalpha^(m-1));

if k==0
    irf = zeros(Nsim,1);
else
    T = 200;
    
    model = arima('Constant',0,'AR',aalpha,'MA',ggamma,'Variance',1);
    
    [Y,~,~] = simulate(model,T,'numPaths',Nsim);

    irf     = zeros(Nsim,1);
    I_k     = [eye(k-1) zeros(k-1,1)];
    
    for i=1:Nsim
        X   = lagY(Y(:,i),k);
        OLS = (X'*X)\(X'*Y(k+1:T,i));
        pphi    = [OLS';I_k];
        irfhold = pphi^m;
        irf(i)= irfhold(1,1);
    end
end

M = mean((irf-ttheta_m).^2);