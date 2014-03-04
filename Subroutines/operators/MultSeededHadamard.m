function [Z] = MultSeededHadamard(X, J, numBlockL, numBlockC, Mblock, Nblock, rp, noBlockError)
% Multiplication of a seeded Hadamard operator with a vector

if (isrow(X) ); X = X'; end
Z = zeros(sum(Mblock), 1);

lastZ = 0;
for l = 1 : numBlockL
    
    Y = 0;
    for c = 1 : numBlockC
        if (J(l, c) ~= 0)
            YY = sqrt(J(l, c) ) * hadamards(X((c - 1) * Nblock + 1 : c * Nblock) );
            Y = Y + YY(rp{l, c}(1 : Mblock(l) ) );
        end
    end
    
    if (max(size(noBlockError) ) > 0); Y(noBlockError{l} ) = -Y(noBlockError{l} ); end
    
    Z(lastZ + 1 : lastZ + Mblock(l) ) = Y;
    lastZ = lastZ + Mblock(l);
end

if (isrow(X) ); Z = Z'; end

end