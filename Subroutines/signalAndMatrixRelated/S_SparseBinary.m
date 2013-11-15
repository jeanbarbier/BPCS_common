function [S] = S_SparseBinary(sze,rho)
% Creates a binary (0,1) signal of density rho of ones

S = [ones(1, floor(sze .* rho)), zeros(1, sze - floor(sze .* rho) )];
S = S(randperm(sze) );

end