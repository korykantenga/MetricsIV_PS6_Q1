function X = lagY(Y,P)
% Create matrix Pith P lags of Y

T = size(Y,1);

W = zeros(T,P);

W(2:T,1) = Y(1:T-1);

for p=1:P-1
    W(p+1:T,p+1) = W(p:T-1,p);
end

X = W(P+1:T,:);