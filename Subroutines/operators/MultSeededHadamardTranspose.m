function [Z] = MultSeededHadamardTranspose(X, J, numBlockL, numBlockC, Mblock, Nblock, rp, noBlockError)
% Multiplication of a seeded Hadamard transposed operator with a vector

if (isrow(X) ); X = X'; end

Z = zeros(numBlockC * Nblock, 1);
ZZ = Z;
zero = zeros(Nblock, 1);

lastZ = 1;
for c = 1 : numBlockC
    
    u = 1; l_ = 1;
    while (J(l_, c) == 0); u = u + Mblock(l_); l_ = l_ + 1; end
    
    Y = 0; YY = 0;
    for (l = 1 : numBlockL)
        
        if (J(l, c) ~= 0)
            XX = zero;
            S = X(u : min(end, u + Nblock - 1) );
            msS = max(size(S) );
            XX(rp{l, c}(1 : min(msS, Mblock(l) ) ) ) = S(1 : min(msS, Mblock(l) ) );
            
            Y = Y + sqrt(J(l, c) ) * hadamards(XX);
            u = u + Mblock(l);
            
            if (max(size(noBlockError) ) > 0)
                XXX = zero;
                XXX(rp{l, c}(noBlockError{l} ) ) = XX(rp{l, c}(noBlockError{l} ) );
                YY = YY + sqrt(J(l, c) ) * hadamards(XXX);
            end
            
        end
    end
    
    Z(lastZ : lastZ + Nblock - 1) = Y;
    
    if (max(size(noBlockError) ) > 0); ZZ(lastZ : lastZ + Nblock - 1) = YY; end
    
    lastZ = lastZ + Nblock;
end

if (max(size(noBlockError) ) > 0); Z = Z - 2 * ZZ; end

if (isrow(X) ); Z = Z'; end

end