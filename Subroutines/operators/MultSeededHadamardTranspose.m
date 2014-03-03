function [Z] = MultSeededHadamardTranspose(X, J, numBlockL, numBlockC, Mblock, Nblock, rp, noBlockError)
% Multiplication of a seeded Hadamard transposed operator with a vector

if (isrow(X) ); X = X'; end

Z = zeros(numBlockC * Nblock, 1);
ZZ = Z;
zero = zeros(Nblock, 1);

lastZ = 1;
for c = 1 : numBlockC
    
    if (c > 2); u = sum(Mblock(1 : c - 2) ) + 1;
    else u = 1; end;
    
    Y = 0; YY = 0;
    l = c - 1;
    while (l <= numBlockL)
        
        if (c == 1); l = l + 1; l = min(l, numBlockL); end
        
        if (J(l, c) == 0); break; end
        
        XX = zero;
        if (numBlockC > 1) S = X(u : min(end, u + Nblock - 1) );
        else S = X(u : end); end
        msS = max(size(S) );
        
        if (numBlockC > 1) XX(rp{l, c}(1 : min(msS, Mblock(l) ) ) ) = S(1 : min(msS, Mblock(l) ) );
        else XX(rp{l, c}(1 : msS ) ) = S(1 : msS); end
        
        Y = Y + sqrt(J(l, c) ) * hadamards(XX);
        u = u + Mblock(l);
        
        if (max(size(noBlockError) ) > 0)
            XXX = zero;
            XXX(rp{l, c}(noBlockError{l} ) ) = XX(rp{l, c}(noBlockError{l} ) );
            YY = YY + sqrt(J(l, c) ) * hadamards(XXX);
        end
        
        l = l + 1;
    end
    
    Z(lastZ : lastZ + Nblock - 1) = Y;
    
    if (max(size(noBlockError) ) > 0) ZZ(lastZ : lastZ + Nblock - 1) = YY; end
    
    lastZ = lastZ + Nblock;
end

if (max(size(noBlockError) ) > 0); Z = Z - 2 * ZZ; end

if (isrow(X) ); Z = Z'; end

end