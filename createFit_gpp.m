function [fitresult, gof] = createFit(a1, b1,cc)
%CREATEFIT(A1,B1)
%  Create a fit.
%
%  Data for 'untitled fit 1' fit:
%      X Input : a1
%      Y Output: b1
%  Output:
%      fitresult : a fit object representing the fit.
%      gof : structure with goodness-of fit info.
%
%  鍙﹁鍙傞槄 FIT, CFIT, SFIT.

%  鐢? MATLAB 浜? 06-Jun-2023 20:43:45 鑷姩鐢熸垚


%% Fit: 'untitled fit 1'.
[xData, yData] = prepareCurveData( a1, b1 );

% Set up fittype and options.
ft = fittype( 'exp(b*(x-c)*(x-c)+d)', 'independent', 'x', 'dependent', 'y' );
opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
opts.Display = 'Off';
b = 25.2;
opts.Lower = [-Inf b -Inf];
opts.StartPoint = [-0.0056 b 2];
opts.Upper = [0 b Inf];
% Fit model to data.
[fitresult, gof] = fit( xData, yData, ft, opts );
if cc==1
% if fitresult.a1 >-9000
% Plot fit with data.
% figure( 'Name', 'untitled fit 1' );
 h = plot( fitresult, xData, yData );
%   h = plot( fitresult);

legend( h, 'Temp vs. ER', 'untitled fit 1', 'Location', 'NorthEast' );
% Label axes
xlabel Temp
ylabel ER
grid on
hold on
% end
end


