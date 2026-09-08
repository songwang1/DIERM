function er = diermER(beta, temperatureC, trefC)
%DIERMER Evaluate ecosystem respiration from the centered DIERM curve.
%   ER = DIERMER(BETA, TEMPERATUREC, TREFC) evaluates
%
%     log(ER) = beta0 + beta1*(Ta-Tref) + beta2*(Ta-Tref)^2.
%
%   BETA may be a 1-by-3 vector or an N-by-3 matrix. For an N-by-3 BETA,
%   TEMPERATUREC must contain N corresponding values. Use
%   diermAnnualMetrics for day-by-location temperature arrays.
%   The linear predictor is clipped to [-50, 50] before exponentiation to
%   prevent numerical overflow. TREFC defaults to 12 degrees C.

if nargin < 3
    trefC = 12;
end
x = double(temperatureC) - trefC;
eta = beta(:,1) + beta(:,2).*x + beta(:,3).*x.^2;
er = exp(max(-50, min(50, eta)));
end
