function y = RMSE(aaalpha,gggamma,Nsim)
% Plot RSME given parameters

k = 0:1:12;
m = 1:1:6;
y = zeros(length(m),length(k));

counter = 0;
for i=1:length(m)
    for j=1:length(k)
        y(i,j) = sqrt(MSE(aaalpha,gggamma,m(i),k(j),Nsim));
        counter = counter + 1 %#ok<NOPRT>
    end
end

end

