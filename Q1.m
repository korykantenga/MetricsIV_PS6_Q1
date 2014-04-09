%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Econometrics IV: Problem Set 6
% Kory Kantenga
% 7 April 2014
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all
clc

profile ON
tic

table1      = 0;
table2      = 1;
rootMSE     = 0;
plotrootMSE = 0;

Nsim   = 1000;
order  = 0:1:12;
aalpha = -0.9:0.2:0.9;
ggamma = -0.9:0.2:0.9;
K      = 12;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%                               TABLE 1
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if table1==1
    tic;
    
    store1  = zeros(10,10,13);
    store2  = zeros(10,10,13);
    
    counter = 0;
    
    for i=1:10
        for j=1:10
            for q=1:13
                if i~=j
                    store1(i,j,q) = ...
                        MSE(aalpha(i),ggamma(j),2,order(q),Nsim);
                    store2(i,j,q) = ...
                        MSE(aalpha(i),ggamma(j),6,order(q),Nsim);
                end
            end
            counter = counter+1 %#ok<NOPTS>
            
        end
    end
    
    [~,table1a] = min(store1,[],3);
    table1a = table1a - ones(size(table1a,1),size(table1a,2));
    
    [~,table1b] = min(store2,[],3);
    table1b = table1b - ones(size(table1b,1),size(table1b,2));
    
    save('table1a','table1a')
    save('table1b','table1b')
    
    toc
    
else
    load('table1b')
    load('table1a')
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%                               RMSE PLOT
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if rootMSE==1
    y1 = RMSE(aalpha(10),ggamma(9),Nsim);
    y2 = RMSE(aalpha(9),ggamma(2),Nsim);
    y3 = RMSE(aalpha(9),ggamma(1),Nsim);
    y4 = RMSE(aalpha(9),ggamma(10),Nsim);
end

if plotrootMSE==1
    figure;
    subplot(2,2,1)
    plot(order,y1)
    title('Root MSE for AR(k) Models \alpha=0.9,\gamma=0.7')
    xlabel('AR(k)')
    ylabel('RMSE')
    legend('m=1','m=2','m=3','m=4','m=5','m=6','Location','Best')
    
    subplot(2,2,2)
    plot(order,y2)
    title('Root MSE for AR(k) Models \alpha=0.7,\gamma=-0.7')
    xlabel('AR(k)')
    ylabel('RMSE')
    legend('m=1','m=2','m=3','m=4','m=5','m=6','Location','Best')
    
    subplot(2,2,3)
    plot(order,y3)
    title('Root MSE for AR(k) Models \alpha=0.7,\gamma=-0.9')
    xlabel('AR(k)')
    ylabel('RMSE')
    legend('m=1','m=2','m=3','m=4','m=5','m=6','Location','Best')
    
    subplot(2,2,4)
    plot(order,y4)
    title('Root MSE for AR(k) Models \alpha=0.7,\gamma=0.9')
    xlabel('AR(k)')
    ylabel('RMSE')
    legend('m=1','m=2','m=3','m=4','m=5','m=6','Location','Best')
    
    print -depsc2 sqrtMSE
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%                               TABLE 2
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Select Model based on Information Criteria

if table2==1
    
    counter = 0;
    
    storeRMSEA = zeros(10,10,2);
    storeRMSEB = zeros(10,10,2);
    
    
    for i=1:length(aalpha)
        for j=1:length(ggamma)
            if i~=j
            [storeRMSEA(i,j,:) storeRMSEB(i,j,:)] = ...
                MSE_IC(aalpha(i),ggamma(j),Nsim,K);
            counter = counter+1 %#ok<NOPTS>
            end
        end
    end
    
end

save('storeRMSEA','storeRMSEA')
save('storeRMSEB','storeRMSEB')

toc
profile OFF
profile VIEWER


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%                               OUTPUT TO LATEX
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%Table 1a: Min MSE (m=2)
s = sym(table1a);
v = vpa(s);
latex(v)

%Table 1b: Min MSE (m=2)
s = sym(table1b);
v = vpa(s);
latex(v)

digits(2);

%Table 2ai: FIC/AIC RMSE (m=2)
s = sym(storeRMSEA(:,:,1));
v = vpa(s);
latex(v)

%Table 2aii: FIC/AIC RMSE (m=6)
s = sym(storeRMSEA(:,:,2));
v = vpa(s);
latex(v)

%Table 2bi: FIC/BIC RMSE (m=2)
s = sym(storeRMSEB(:,:,1));
v = vpa(s);
latex(v)

%Table 2bii: FIC/BIC RMSE (m=6)
s = sym(storeRMSEB(:,:,2));
v = vpa(s);
latex(v)


