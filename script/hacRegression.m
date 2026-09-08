function fit = hacRegression(y,X,maxLag)
%HACREGRESSION Fit OLS with Newey-West heteroskedasticity/autocorrelation SE.
if nargin<3, maxLag=3; end
y=double(y(:)); X=double(X); ok=isfinite(y)&all(isfinite(X),2); y=y(ok); X=X(ok,:);
n=size(X,1); k=size(X,2); beta=X\y; residual=y-X*beta;
S=zeros(k);
for t=1:n, S=S+(residual(t)^2).*(X(t,:)'*X(t,:)); end
for lag=1:maxLag
    weight=1-lag/(maxLag+1); gamma=zeros(k);
    for t=lag+1:n
        gamma=gamma+residual(t)*residual(t-lag).*(X(t,:)'*X(t-lag,:));
    end
    S=S+weight.*(gamma+gamma');
end
bread=pinv(X'*X); covariance=bread*S*bread;
se=sqrt(max(0,diag(covariance))); statistic=beta./se;
p=2*(1-tcdf(abs(statistic),max(1,n-k)));
fit=struct('beta',beta,'se',se,'p',p,'covariance',covariance,'residual',residual,'n',n);
end
