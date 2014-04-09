function irf = arIRF(m,k,beta)
% Recursively calculate impulse response for Gaussian AR(p) process
if k==0
    irf = zeros(k,1);
else
    Y      = zeros(m+k+1,1);
    Y(k+1) = 1; % one sd shock
    
    for i=k+2:k+m+1
        Y(i) = Y(i-k:i-1)'*beta;
    end
    
    irf = Y(m+1:m+k);
end

end

