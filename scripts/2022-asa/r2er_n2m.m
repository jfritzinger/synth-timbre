%% Explainable variance explained 
function [hat_r2er, r2] = r2er_n2m(x, y)
%  Neuron to model approx. unbiased estimator of r2er.
%    Parameters
%    ----------
%    x : numpy.ndarray
%    m observations model predictions
%    y : numpy.ndarray
%    N neurons X n reps X m observations of data
%    
%    Returns
%    -------
%    r2er : an estimate of the r2 between model and expected value of data
%    --------

[n, m] = size(y);
% estimate of trial-to-trial variability
sig2_hat = mean(var(y));
% mean center the model
x_ms = x - mean(x);
% get average across reps (rows)
y = mean(y);
% mean center data average
y_ms = y - mean(y);

% dot product of mean centered model and data
% are biased numerator
xy2 = sum(x_ms.*y_ms).^2;

% individual variances of model and data (columns)
x2 = sum(x_ms.^2);
y2 = sum(y_ms.^2);

% biased denominator
x2y2 = x2*y2;
r2 = xy2/x2y2;

% unbias numerator and denominator
ub_xy2 = xy2 - sig2_hat/n * x2;
ub_x2y2 = x2y2 - (m-1)*sig2_hat/n*x2;

% form ratio of individually unbiased estimates from r^2_er
hat_r2er = ub_xy2/ub_x2y2;
end