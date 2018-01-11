function [corr] =  autocor(X, maxlag)
N = length(X); 
corr = zeros(maxlag,1); 
mu = mean(X);
var = sum((X-mu*ones(N,1)).^2)/(N-1); 
for t = 1:maxlag
    sums = 0; 
    for s = 1:N-t
        term = (X(s)-mu)*(X(s+t)-mu);
        sums = sums + term; 
    end
    corr(t) = sums/((N-t)*var);
end