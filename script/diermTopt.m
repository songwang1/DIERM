function toptC = diermTopt(beta, trefC)
%DIERMTOPT Return the vertex of the fitted log-ER parabola.
%   Predicted non-negative beta2 values are forced to -1e-10 before the
%   vertex is calculated, preserving the unimodal response assumption.

if nargin < 2
    trefC = 12;
end
beta2 = min(double(beta(:,3)), -1e-10);
toptC = trefC - double(beta(:,2))./(2.*beta2);
end
