function [G] = createSeededHadamardMat(Mblock, Nblock, numBlockL, numBlockC, J, rp, noBlockError)
% creates a block seeded matrix from the variance matrix J

G = zeros(sum(Mblock), numBlockC * Nblock);
h = hadamard(Nblock);
startL = 1; stopL = Mblock(1);

for l = 1 : numBlockL
    
    for c = 1 : numBlockC
        if (max(size(rp{l, c} ) ) > 0)
            hh = sqrt(J(l, c) ) .* h(rp{l, c}, :);
            G(startL : stopL, Nblock .* (c - 1) + 1 : Nblock .* c) = hh(1 : Mblock(l), :);
            
            if (max(size(noBlockError) ) > 0);
                GG = G(startL : stopL, Nblock .* (c - 1) + 1 : Nblock .* c);
                GG(noBlockError{l}, :) = -GG(noBlockError{l}, :);
                G(startL : stopL, Nblock .* (c - 1) + 1 : Nblock .* c) = GG;
            end
            
        else
            G(startL : stopL, Nblock .* (c - 1) + 1 : Nblock .* c) = 0;
        end
    end
    
    if (l < numBlockL); startL = stopL + 1; stopL = stopL + Mblock(l + 1); end
    if (l == numBlockL); startL = stopL + 1; stopL = sum(Mblock); end
end

end