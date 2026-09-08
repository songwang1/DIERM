function fit = clusterRobustRegression(y,X,cluster,weights)
%CLUSTERROBUSTREGRESSION Fit WLS with cluster-robust sandwich covariance.
if nargin<4, weights=ones(size(y)); end
y=double(y(:)); X=double(X); cluster=cluster(:); weights=double(weights(:));
ok=isfinite(y)&all(isfinite(X),2)&isfinite(weights)&weights>0;
y=y(ok); X=X(ok,:); cluster=cluster(ok); weights=weights(ok);
root=sqrt(weights); Xw=X.*root; yw=y.*root; beta=Xw\yw;
residual=(y-X*beta).*root; bread=pinv(Xw'*Xw); meat=zeros(size(X,2));
groups=unique(cluster);
for g=groups'
    use=cluster==g; score=Xw(use,:)'*residual(use); meat=meat+score*score';
end
n=numel(y); k=size(X,2); G=numel(groups);
correction=(G/(G-1))*((n-1)/(n-k)); covariance=correction*bread*meat*bread;
se=sqrt(max(0,diag(covariance))); fit=struct('beta',beta,'se',se,'covariance',covariance,'n',n,'clusters',G);
end
