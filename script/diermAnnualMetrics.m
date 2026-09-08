function out = diermAnnualMetrics(temperatureC, beta, trefC)
%DIERMANNUALMETRICS Calculate Topt, HOT days, and annual ER.
%   TEMPERATUREC is day-by-location. BETA is location-by-3. Annual ER is
%   returned in g C m^-2 yr^-1, assuming daily ER is in micromol CO2
%   m^-2 s^-1.

if nargin < 3
    trefC = 12;
end
conversion = 86400 * 12e-6;
topt = diermTopt(beta, trefC);
hot = sum(temperatureC > reshape(topt, 1, []), 1, 'omitnan')';
x = double(temperatureC) - trefC;
eta = reshape(beta(:,1),1,[]) + reshape(beta(:,2),1,[]).*x + ...
    reshape(min(beta(:,3),-1e-10),1,[]).*x.^2;
dailyER = exp(max(-50, min(50, eta)));
annualER = sum(dailyER,1,'omitnan')'.*conversion;
invalid = ~all(isfinite(beta),2) | ~isfinite(topt);
hot(invalid) = NaN;
annualER(invalid) = NaN;
out = table(topt, hot, annualER, ...
    'VariableNames', {'topt_c','hot_days','annual_er_g_c_m2'});
end
