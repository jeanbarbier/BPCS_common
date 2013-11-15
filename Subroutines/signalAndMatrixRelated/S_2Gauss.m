function [S] = S_2Gauss(size,rho,m_small,m_big,var_small,var_big)
% Creates an homogeneous random Bi-Gaussian vector with rho as fraction
% of big components, var_big as expectation for this gaussian and
% var_big for it's variance. var_small is the variance of the small
% components of the produced signal and m_small its average.
% size is the number of components of the original signal

num_non_zero = floor(rho .* size);
num_zero = size - num_non_zero;

SS = ([randn(1, num_zero) .* sqrt(var_small) + m_small, randn(1, num_non_zero) .* sqrt(var_big) + m_big] );
S(randperm(size) ) = SS;

end