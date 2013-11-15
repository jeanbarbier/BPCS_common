function [S] = S_SparseGaussComplex(N, rho, mGauss, varGauss)
% Creates a complex Gauss-Bernoulli complex signal of distribution p(x) ~ (1 - rho) * delta(x) + rho / sqrt(2 * pi * var_gauss) * exp(-(|x| - m_gauss)^2 / (2 * var_gauss) )

SnonZero = mGauss + (randn(1, floor(rho .* N) ) + 1i .* randn(1, floor(rho .* N) ) ) .* sqrt(varGauss);
S = zeros(1, N);
S(randperm(N, floor(rho .* N) ) ) = SnonZero;

end