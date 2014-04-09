function [AIC BIC] = AICBIC(k,data)

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

end

