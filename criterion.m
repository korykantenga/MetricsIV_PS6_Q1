function [AIC BIC FIC] = criterion(m,k,K,data)

% Output AIC, BIC and FIC for lag k, horizon m and data vector (Tx1)
% K is the size of the full model for FIC

T  = size(data,1);

X   = lagY(data,k);
OLS = (X'*X)\(X'*data(k+1:T));
if k>0
    RSS = (data(k+1:T)-X*OLS)'*(data(k+1:T)-X*OLS);
else
    RSS = data(:)'*data(:);
end

AIC = -1*RSS - 2*k;
BIC = -1*RSS - k*log(T-k);

XFull  = lagY(data,K);
BFull  = (XFull'*XFull)\(XFull'*data(K+1:T));
SSigma = ((data(K+1:T)-XFull*BFull)'*(data(K+1:T)-XFull*BFull))/...
    (T-2*K);

Sm     = zeros(K,k);
Sm(1:k,1:k) = eye(k);

if k>0
    IRF    = arIRF(m,k,OLS);
    D      = [IRF;zeros(K-k,1)];
else
    IRF    = arIRF(m,k,OLS);
    D = [IRF;zeros(K-1,1)];
end

if k>0
    FIC = (D'*(BFull-Sm*OLS))^2+2*(SSigma*D'*Sm)*((X'*X)\Sm'*D);
else
    FIC = 0;
end


end

