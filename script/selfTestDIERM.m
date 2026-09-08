function selfTestDIERM
%SELFTESTDIERM Run lightweight tests without external data files.
rng(7);
temperature=linspace(-5,35,30)'; trueBeta=[1.2 0.08 -0.006];
respiration=diermER(trueBeta,temperature,12).*(1+0.03*randn(size(temperature)));
fit=fitSiteYearBeta(temperature,respiration,12);
assert(~isempty(fit)&&fit.beta2<0,'The constrained site-year fit failed.');
assert(abs(fit.implied_topt_c-diermTopt(trueBeta,12))<3,'Recovered Topt is unexpectedly inaccurate.');

X=randn(80,5); Y=[X(:,1)+0.2*X(:,2),X(:,3).^2,-0.01-0.002*abs(X(:,4))];
model=fitExtraTreesRegressor(X,Y,NumTrees=20,MinLeafSize=4,MaxFeatures=0.8,Seed=11);
prediction=predictExtraTreesRegressor(model,X);
assert(isequal(size(prediction),size(Y))&&all(isfinite(prediction),'all'),'Extra-Trees prediction failed.');

dailyTemperature=repmat(linspace(-5,35,365)',1,2);
annual=diermAnnualMetrics(dailyTemperature,[trueBeta;trueBeta],12);
assert(all(annual.hot_days>=0&annual.hot_days<=365),'HOT-day counts are invalid.');
fprintf('DIERM MATLAB self-test passed.\n');
end
