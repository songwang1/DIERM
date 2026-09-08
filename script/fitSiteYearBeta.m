function fit = fitSiteYearBeta(temperatureC, respiration, trefC)
%FITSITEYEARBETA Fit one constrained DIERM curve to temperature-bin means.
%   The function uses multiple starting values and bounded nonlinear least
%   squares. beta2 is constrained to [-2, -1e-10], which imposes a
%   downward-opening log-response curve.

if nargin < 3
    trefC = 12;
end
t = double(temperatureC(:));
y = double(respiration(:));
ok = isfinite(t) & isfinite(y) & y > 0;
t = t(ok); y = y(ok);
if numel(y) < 3
    fit = table();
    return
end
x = t - trefC;
scale = max([std(y), mean(y)*0.1, 1e-6]);
p = polyfit(x, log(max(y,1e-8)), 2);
q0 = min(19.9,max(-19.9,p(3)));
q1 = min(4.9,max(-4.9,p(2)));
q2 = min(-1e-8,max(-1.9,p(1)));
curvatures = [-0.05 -0.01 -0.002 -0.0002 -1e-7];
starts = [q0 q1 q2];
for c = curvatures
    starts(end+1,:) = [log(max(mean(y),1e-8)), q1, c]; %#ok<AGROW>
    starts(end+1,:) = [log(max(max(y),1e-8)), 0, c]; %#ok<AGROW>
end
lb = [-20 -5 -2]; ub = [20 5 -1e-10];
objective = @(b) (exp(max(-50,min(50,b(1)+b(2).*x+b(3).*x.^2)))-y)./scale;
options = optimoptions('lsqnonlin','Display','off','MaxFunctionEvaluations',5000);
bestSSE = inf; bestBeta = nan(1,3); success = false;
for k = 1:size(starts,1)
    start = max(lb+1e-10,min(ub-1e-10,starts(k,:)));
    try
        [candidate,~,~,exitflag] = lsqnonlin(objective,start,lb,ub,options);
        prediction = exp(max(-50,min(50,candidate(1)+candidate(2).*x+candidate(3).*x.^2)));
        sse = sum((prediction-y).^2);
        if sse < bestSSE
            bestSSE = sse; bestBeta = candidate; success = exitflag > 0;
        end
    catch
        % Continue to the next starting point if one optimization fails.
    end
end
if ~isfinite(bestSSE)
    fit = table();
    return
end
sst = sum((y-mean(y)).^2);
r2 = 1-bestSSE/sst;
if sst == 0, r2 = NaN; end
if bestBeta(3) < -1e-8
    impliedTopt = trefC-bestBeta(2)/(2*bestBeta(3));
    impliedC = bestBeta(1)-bestBeta(2)^2/(4*bestBeta(3));
else
    impliedTopt = NaN; impliedC = NaN;
end
fit = table(bestBeta(1),bestBeta(2),bestBeta(3),r2, ...
    sqrt(bestSSE/numel(y)),impliedTopt,impliedC,abs(bestBeta(3))<1e-5,success, ...
    'VariableNames',{'beta0','beta1','beta2','r2','rmse_binned', ...
    'implied_topt_c','implied_c','beta2_near_zero','optimizer_success'});
end
