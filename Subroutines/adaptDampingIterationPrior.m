% adaptive damping procedure

opt.dump_mes = 0.;

if (t==1)
    priorTry.R = var_2 ./ var_1 + prior.av_mess;
    priorTry.S2 = 1 ./ var_1;
    
    % Computing mean and variance with a given Prior
    priorTry = F(priorTry);
    
    % free entropy
    if (opt.save_speed ~= 1); loglike = -0.5 * measure_rate * (mean2(((Y - (G * priorTry.av_mess.').').^2) ./ n_and_e.var_noise + log(n_and_e.var_noise + (G.^2 * priorTry.var_mess.').') + log(2 * pi) ) );
    else loglike = -0.5 * measure_rate * (mean2(((Y - (G * priorTry.av_mess.').').^2) ./ n_and_e.var_noise + log(n_and_e.var_noise + (G2 * priorTry.var_mess.').') + log(2 * pi) ) ); end
    
    % Kulback Leibler divergence
    Part_1 = (priorTry.var_mess + (priorTry.av_mess - priorTry.R).^2) ./ (2 .* priorTry.S2);
    m = 0; s2 = 1;
    Z_tilde = (1 - prior.rho) * exp(-priorTry.R.^2 ./ (2 * priorTry.S2) ) + prior.rho * sqrt(priorTry.S2 ./ (priorTry.S2 + s2) ) .* exp(-((priorTry.R - m).^2) ./ (2 * (priorTry.S2 + s2) ) );
    Part_2 = log(Z_tilde);
    DL = Part_1 + Part_2;
    FreeEntropy = loglike  + sum(DL) / N;
    
else
    
    notok=1;
    while (notok == 1)
        priorTry.R = (1 - opt.dump_mes) * (var_2 ./ var_1 + prior.av_mess) + opt.dump_mes * prior.R;
        priorTry.S2 = (1 - opt.dump_mes) * (1 ./ var_1) + opt.dump_mes * prior.S2;
        
        %Computing mean and variance with a given Prior
        priorTry = F(priorTry);
        
        % free entropy
        if (opt.save_speed ~= 1); loglike = -0.5 * measure_rate * (mean2(((Y - (G * priorTry.av_mess.').').^2) ./ n_and_e.var_noise + log(n_and_e.var_noise + (G.^2 * priorTry.var_mess.').') + log(2 * pi) ) );
        else loglike = -0.5 * measure_rate * (mean2(((Y - (G * priorTry.av_mess.').').^2) ./ n_and_e.var_noise + log(n_and_e.var_noise + (G2 * priorTry.var_mess.').') + log(2 * pi) ) ); end
        
        % Kulback Leibler divergence
        Part_1 = (priorTry.var_mess + (priorTry.av_mess - priorTry.R).^2) ./ (2 .* priorTry.S2);
        m = prior.param_1; s2 = prior.param_2;
        Z_tilde = (1 - prior.rho) * exp(-priorTry.R.^2 ./ (2 * priorTry.S2) ) + prior.rho * sqrt(priorTry.S2 ./ (priorTry.S2 + s2) ) .* exp(-((priorTry.R - m).^2) ./ (2 * (priorTry.S2 + s2) ) );
        Part_2 = log(Z_tilde);
        DL = Part_1 + Part_2;
        FreeEntropy = loglike  + sum(DL) / N;
        
        if ((FreeEntropy - FreeEntropy_OLD) / abs(FreeEntropy_OLD) > -0.2); notok=0;
        else opt.dump_mes = opt.dump_mes + (1 - opt.dump_mes) / 2; end
        
        if (opt.dump_mes > 1 - 1e-6); notok=0; end
    end
end

FreeEntropy_OLD = FreeEntropy;
prior = priorTry;
V = V_new; W = W_new;

% free entropy
if (opt.save_speed ~= 1); loglike = -0.5 * measure_rate * (mean2(((Y - (G * prior.av_mess.').').^2) ./ n_and_e.var_noise + log(n_and_e.var_noise + (G.^2 * prior.var_mess.').') + log(2 * pi) ) );
else loglike = -0.5 * measure_rate * (mean2(((Y - (G * prior.av_mess.').').^2) ./ n_and_e.var_noise + log(n_and_e.var_noise + (G2 * prior.var_mess.').') + log(2 * pi) ) ); end

% Kulback Leibler divergence
Part_1 = (prior.var_mess + (prior.av_mess - prior.R).^2) ./ (2 .* prior.S2);
m = prior.param_1; s2 = prior.param_2;
Z_tilde = (1 - prior.rho) * exp(-prior.R.^2 ./ (2 * prior.S2) ) + prior.rho * sqrt(prior.S2 ./ (prior.S2 + s2) ) .* exp(-((prior.R - m).^2) ./ (2 * (prior.S2 + s2) ) );
Part_2 = log(Z_tilde);
DL = Part_1 + Part_2;
FreeEntropy = loglike  + sum(DL) / N;