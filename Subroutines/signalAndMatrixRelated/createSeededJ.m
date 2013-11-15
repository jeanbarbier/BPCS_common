function [J] = createSeededJ(numBlockL, numBlockC, JJ, w, N)
% creates a seeded variance matrix with one upper diagonal block of variance JJ / N and w sub-diagonal blocks of variance 1 / N

J = zeros(numBlockL, numBlockC);
ww = [0 : w];

for l = 1 : numBlockL
    
    for c = 1 : numBlockC
        if (c == l + 1);
            J(l, c) = JJ ./ N;
        elseif (sum(c == l - ww) > 0)
            J(l, c) = 1 ./ N;
        end
        
    end
    
end

end