function [Z] = MultSeededHadamardSquarred(X, J, numBlockL, numBlockC, Mblock, Nblock)
% Multiplication of a squarred seeded Hadamard operator with a vector

if (isrow(X) ); X = X'; end
Z = zeros(sum(Mblock), 1);

lastZ = 0;
for l = 1 : numBlockL
    
    Y = 0;
    for c = 1 : numBlockC
        if (J(l, c) ~= 0)
            YY = J(l, c) * sum(X((c - 1) * Nblock + 1 : c * Nblock) ) * ones(Mblock(1), 1);
            Y = Y + YY(1 : Mblock(l) );
        end
    end
    
    Z(lastZ + 1 : lastZ + Mblock(l) ) = Y;
    lastZ = lastZ + Mblock(l);
end

if (isrow(X) ); Z = Z'; end

end