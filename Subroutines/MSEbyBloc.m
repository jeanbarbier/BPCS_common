function [MSE] = MSEbyBloc(X, S, numBlockC, Nblock)
% compute the MSE by block for the seeded matrices

MSE = zeros(1, numBlockC);

if (isreal(S) && (numBlockC > 1) )
    for l = 1 : numBlockC; MSE(l) = mean((X((l - 1) * Nblock + 1: l * Nblock) - S((l - 1) * Nblock + 1: l * Nblock) ).^2); end
elseif ((isreal(S) == 0) && (numBlockC > 1) )
    for l = 1 : numBlockC; MSE(l) = 0.5 .* (mean((real(X((l - 1) * Nblock + 1: l * Nblock) ) - real(S((l - 1) * Nblock + 1: l * Nblock) ) ).^2) + mean((imag(X((l - 1) * Nblock + 1: l * Nblock) ) - imag(S((l - 1) * Nblock + 1: l * Nblock) ) ).^2) ); end
end

end