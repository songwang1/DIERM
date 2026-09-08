function out = weightedSummary(values,weights,prefix)
%WEIGHTEDSUMMARY Return an area-weighted mean and weighted quantiles.
v=double(values(:)); w=double(weights(:));
ok=isfinite(v)&isfinite(w)&w>0; v=v(ok); w=w(ok);
out=struct();
if isempty(v)
    q=[NaN NaN NaN]; mu=NaN;
elseif isscalar(v)
    q=repmat(v,1,3); mu=v;
else
    [v,order]=sort(v); w=w(order); mu=sum(v.*w)/sum(w);
    cdf=(cumsum(w)-0.5*w)/sum(w);
    [cdf,keep]=unique(cdf); v=v(keep);
    q=interp1(cdf,v,[0.05 0.5 0.95],'linear','extrap');
end
out.(prefix+"_mean")=mu; out.(prefix+"_p05")=q(1);
out.(prefix+"_median")=q(2); out.(prefix+"_p95")=q(3);
end
