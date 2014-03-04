function [Z] = MultSeededHadamardTransposeSquarred(X, J, numBlockL, numBlockC, Mblock, Nblock)
% Multiplication of a squarred seeded Hadamard transposed operator with a vector

if (isrow(X) ); X = X'; end
Z = zeros(numBlockC * Nblock, 1);

lastZ = 1;
for c = 1 : numBlockC
    
    u = 1; l_ = 1;
    while (J(l_, c) == 0); u = u + Mblock(l_); l_ = l_ + 1; end
    
    Y = 0;
    for (l = 1 : numBlockL)
        
        if (J(l, c) ~= 0)
            Y = Y + J(l, c) * sum(X(u : min(end, u + Mblock(l) - 1) ) ) * ones(Nblock, 1);
            u = u + Mblock(l);
        end
    end
    
    Z(lastZ : lastZ + Nblock - 1) = Y;
    lastZ = lastZ + Nblock;
end

if (isrow(X) ); Z = Z'; end

end