function [G] = createSeededRandomMatrix(type, J, Mblock, Nblock)
% Creates a random seeded real gaussian or seeded complex matrix from the
% J variance matrix.

[numBlockL, numBlockC] = size(J);
N = Nblock .* numBlockC;

if (strcmp(type, 'real') ); G = randn(sum(Mblock), N);
else G = exp(1i .* 2 .* pi .* rand(sum(Mblock), N) ); end

% "seeding" of the matrix
startL = 1; stopL = Mblock(1);
for l = 1 : numBlockL
    for c = 1 : numBlockC; G(startL : stopL, Nblock .* (c - 1) + 1 : Nblock .* c) = G(startL : stopL, Nblock .* (c - 1) + 1 : Nblock .* c) .* sqrt(J(l, c) ); end
    if (l < numBlockL); startL = stopL + 1; stopL = stopL + Mblock(l + 1); end
    if (l == numBlockL); startL = stopL + 1; stopL = sum(Mblock); end
end

end