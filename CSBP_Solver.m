function [X, varargout] = CSBP_Solver(Y, G, opt)
check_size();
set_parameters();
display_information();
initialisation();

% Construction of the prior dependent class and noise_and_error
prior = Prior(opt.signal_rho, N, measure_rate, opt.learn, opt.prior, opt.dump_learn, R_init, S2_init, av_mess_init, var_mess_init, opt.method, prior_param{1}, prior_param{2}, prior_param{3}, prior_param{4}); F = str2func(prior.func);
n_and_e = noise_and_error(opt.conv, opt.var_noise, opt.dump_learn);

t = 1; print_to_screen();
% initial error (should be rho)
if (max(size(opt.signal) ) > 2); n_and_e = n_and_e.compute_true_MSE(opt.signal,prior.av_mess); MSEt(t) = n_and_e.true_error; end
% initial MSE by block
if ((opt.MSEbyBlock > 0) && (mod(t, opt.MSEbyBlock) == 0) && (opt.numBlockC > 1) ); MSEblocks{t} = MSEbyBloc(prior.av_mess, opt.signal, opt.numBlockC, opt.Nblock); end

figure;
% Starting main code
while (t <= opt.nb_iter)
    
    switch (opt.method)
        case ('AMP'); AMP();
        case ('AMPcomplex'); AMPcomplex();
        case ('AMPtap'); AMPtap();
        case ('AMPseededHadamard'); AMPseededHadamard();
        case ('AMPseededFourier'); AMPseededFourier();
    end
    
    % Test of the convergence
    n_and_e = n_and_e.compute_convergence(prior.av_mess_old,prior.av_mess);
    if ((n_and_e.convergence < opt.conv) && (t > 10) ); fprintf('Converged, convergence = %e',n_and_e.convergence); break; end;
    
    % Test of reconstruction on the fly knowing the original signal
    if (max(size(opt.signal) ) > 2)
        n_and_e = n_and_e.compute_true_MSE(opt.signal,prior.av_mess);
        if ((n_and_e.true_error < opt.conv) ); fprintf('Solution found, true error = %e',n_and_e.true_error); break; end;
    end
    
    % Learning of the noise if activated
    if (opt.option_noise == 1); n_and_e = n_and_e.learn_noise(Y,W_new,V_new); end
    
    % print infos to screen
    if ((opt.print > 0) && (mod(t, opt.print) == 0) ); print_to_screen(); end
    
    % MSE by block
    if ((opt.MSEbyBlock > 0) && (mod(t, opt.MSEbyBlock) == 0) && (opt.numBlockC > 1) )
        MSEblocks{t + 1} = MSEbyBloc(prior.av_mess, opt.signal, opt.numBlockC, opt.Nblock);
        semilogy(MSEblocks{t + 1} );
        if (t == opt.MSEbyBlock); MSEmax = max(max(MSEblocks{t + 1} ) ); end
        axis ([1, opt.numBlockC, 0, MSEmax] ); drawnow;
    end
    
    % MSE as a function of the iterations, to be compared with density evolution
    if (max(size(opt.signal) ) > 2)
        MSEt(t + 1) = n_and_e.true_error;
        if ((isnan(n_and_e.true_error) == 1) || (n_and_e.true_error > 1e5) ); error('The algorithm did not converge'); end
    end
    
    t = t + 1;
    
end

close(gcf);

% can output the MSE a function of iterations
varargout{1} = prior;
varargout{2} = n_and_e;
if (max(size(opt.signal) ) > 2); varargout{3} = MSEt; end
if ((opt.MSEbyBlock > 0) && (opt.numBlockC > 1) ); varargout{4} = MSEblocks; end

X = prior.av_mess;

end