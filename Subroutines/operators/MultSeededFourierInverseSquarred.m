function [Z] = MultSeededFourierInverseSquarred(X, J, numBlockL, numBlockC, Mblock, Nblock)
% Multiplication of a squarred seeded Fourier operator transposed with a vector

if (isrow(X) ); X = X.'; X = [X; zeros(2 * Nblock, 1)];
else X = [X; zeros(2 * Nblock, 1)]; end
Z = zeros(numBlockC * Nblock, 1);

lastZ = 1;
for c = 1 : numBlockC
    
    u = 1;
    Y = 0;
    for l = 1 : numBlockL
        
        if (J(l, c) ~= 0); Y = Y + J(l, c) * sum(X(u : u + Mblock(l) - 1) ) * ones(Nblock, 1); end;
        
        u = u + Mblock(l);
    end
    
    Z(lastZ : lastZ + Nblock - 1) = Y;
    lastZ = lastZ + Nblock;
end

if (isrow(X) ); Z = Z.'; end

end