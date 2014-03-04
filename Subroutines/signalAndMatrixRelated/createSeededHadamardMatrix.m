function [G] = createSeededHadamardMatrix(J, Mblock, Nblock, rp, flipedSigns)
% Creates a full seeded Hadamard J variance matrix.

[numBlockL, numBlockC] = size(J);
N = Nblock .* numBlockC;

G = zeros(sum(Mblock), N);
H = hadamard(Nblock);

% "seeding" of the matrix
startL = 1; stopL = Mblock(1);
for l = 1 : numBlockL
    for c = 1 : numBlockC
        if (J(l, c) ~= 0);
            HH = H;
            HH = HH(rp{l, c}(1 : Mblock(l) ), :);
            G(startL : stopL, Nblock * (c - 1) + 1 : Nblock * c) = HH .* sqrt(J(l, c) );
        end
        
        if (max(size(flipedSigns) ) > 0)
            GG = G(startL : stopL, Nblock * (c - 1) + 1 : Nblock * c);
            GG(flipedSigns{l}, :) = -GG(flipedSigns{l}, :);
            G(startL : stopL, Nblock * (c - 1) + 1 : Nblock * c) = GG;
        end
    end
    
    if (l < numBlockL); startL = stopL + 1; stopL = stopL + Mblock(l + 1); end
    if (l == numBlockL); startL = stopL + 1; stopL = sum(Mblock); end
end

end