function [rperm, flipedSigns] = createRandomLinesAndSignsPermutationForOperators(numBlockC, numBlockL, J, Mblock, Nblock)
% creates the random permutations of the lines for the fourier or hadamard operators avoiding repetitions and the first (only ones) mode

% selection of the lines of the matrices that are sign-switched
for (l = 1 : numBlockL); flipedSigns{l} = randperm(Mblock(l), floor(Mblock(l) ./ 2) ); end

% selection of the operator modes for each column of the matrices
for (c = 1 : numBlockC)
    
    fullPerm = randperm(Nblock);
    % No first mode
    fullPerm(fullPerm == 1) = fullPerm(end);
    
    rpC = fullPerm(1 : end - 1);
    
    u = 1;
    for (l = 1 : numBlockL)
        if (J(l, c) ~= 0)
            
            if (u >= max(size(rpC) ) - Mblock(1) );
                u = 1;
                rpC = rpC(randperm(Nblock - 1) );
            end
            
            rperm{l, c} = rpC(u : u + Mblock(l) - 1);
            u = u + Mblock(l);
        else
            rperm{l, c} = [];
        end
    end
    
end

end