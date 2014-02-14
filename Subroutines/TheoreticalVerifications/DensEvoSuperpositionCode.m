function [E_new, varargout] = DensEvoSuperpositionCode(B, alpha, Delta, E, sizeMC, type)
% Density evolution iteration for the superposition codes on the Nishimori line for a scaling F ~ O(1 / sqrt(L) )

minn = -10; maxx = -minn; accuracyInt = 1e-10;

if (B == 2)
    
    S2 = (Delta / B + E) ./ alpha;
    Rnd = randn(B, sizeMC) ./ sqrt(S2);
    
    E_new = integral(@(z) exp(-0.5 .* z.^2 ) ./ sqrt(2 .* pi) .* (1 + exp(1 ./ S2 + z .* sqrt(2 ./ S2) ) ).^(-2), minn, maxx, 'AbsTol', accuracyInt);
    
    varargout{1} = 0;
    
    if (strcmp(type, 'sectionErrorRate') )
        % section error rate
        fa_1 = (1 + exp(-1 ./ S2 - Rnd(1, :) + Rnd(2, :) ) ).^(-1);
        cond1 = fa_1 < 0.5;
        
        varargout{1} = mean(cond1);
    end
    
elseif (B > 2)
    
    % B integrals by MC...
    S2 = (Delta / B + E) ./ alpha;
    Rnd = randn(B, sizeMC) ./ sqrt(S2);
    
    % MSE
    intA = mean((1 + (exp(-1 ./ S2 - Rnd(1, :) ) .* sum(exp(Rnd(2 : end, :) ) ) ).^(-1) ).^(-2) );
    intC = mean((1 + exp(1 ./ S2 + Rnd(1, :) - Rnd(2, :) ) + exp(-Rnd(2, :) ) .* sum(exp(Rnd(3 : end, :) ) ) ).^(-2) );
    
    varargout{1} = 0;
    
    if (strcmp(type, 'sectionErrorRate') )
        
        % section error rate
        fa_1 = (1 + exp(-1 ./ S2 - Rnd(1, :) ) .* sum(exp(Rnd(2 : end, :) ) ) ).^(-1);
        
        for bb = 2 : B;
            fa_0(bb - 1, 1 : sizeMC) = (1 + exp(1 ./ S2 + Rnd(1, :) - Rnd(bb, :) ) + exp(-Rnd(bb, :) ) .* sum(exp(Rnd(setxor([2 : B], bb), :) ) ) ).^(-1);
            cond1(bb - 1, 1 : sizeMC) = fa_0(bb - 1, :) > fa_1;
        end
        
        cond2 = mean(cond1) > 0;
        
        secErrRate = mean(cond2);
        varargout{1} = secErrRate;
    end
    
    E_new = 1 ./ B .* intA + (1 - 1 ./ B) .* intC;
    
end

end